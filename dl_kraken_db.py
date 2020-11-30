#!/usr/bin/env python

"""
Downloads data from ensembl genomes for creawting a kraken database, allowing considerably greater
coverage of fungal species than the standard kraken database.

Run `dl_kraken_db.py --help` for usage information
"""

import argparse
import ftplib
import glob
import gzip
import ntpath
import os.path
import pickle
import pprint
import re
import subprocess
import shutil
import time
import sys
import io
from bs4 import BeautifulSoup
from pathlib import Path
from ftplib import FTP
from joblib import Parallel, delayed
from Bio import SeqIO

def main():
    
    global data_dir, metadata_files, retryCount, retrySleep

    # FTP download retry parameters
    retryCount=5
    retrySleep=5

    metadata_files={
        'bacteria': 'species_metadata_EnsemblBacteria',
        'fungi': 'species_metadata_EnsemblFungi'
    }

    parser = argparse.ArgumentParser(
        description="Create ensembl-based kraken database.", 
        epilog="In the event of download failures, the script can be rerun \
            and will continue to download missing records without \
            redownloading successfuly retrieved files.",
        add_help=True
    )
    parser.add_argument('--db', help='Path to kraken database',required=True)
    parser.add_argument('--cleanup', help='Deletes downloaded data and exits',type=bool)
    args = parser.parse_args()
    data_dir = args.db
    
    if args.cleanup:
        cleanup()

    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    if not os.path.exists('{}/tmp'.format(data_dir)):
        os.mkdir('{}/tmp'.format(data_dir))

    if not os.path.exists('{}/ensembl'.format(data_dir)):
        os.mkdir('{}/ensembl'.format(data_dir))

    if not os.path.exists('{}/flags'.format(data_dir)):
        os.mkdir('{}/flags'.format(data_dir))

    for div in ('fungi','bacteria'):
        print('Processing {}...'.format(div))
        metadata_xml = metadata_files[div] + '.xml'
        metadata_pkl = metadata_files[div] + '.pkl'
        
        # parsing the metadata files takes a while, so outputs are pickled to save time on reruns...
        if os.path.exists('{}/{}'.format(data_dir,metadata_pkl)):
            print('Unpickling metadata...')
            genomes=pickle.load( open( "{}/{}".format(data_dir,metadata_pkl), "rb" ) )
        else:
            download_metadata(metadata_xml,div)
            genomes = parse_metadata(metadata_xml)

        print('Downloading genomes...')
        Parallel(n_jobs=4, verbose=5)(delayed(
            download_genome)(genome,genomes[genome]['name'],genomes[genome]['taxid'],div)
        for genome in genomes.keys()) 


def parse_metadata(metadata_file):

    """
    Parses metadata xml files downloaded from ensembl.

    The ensembl ftp site organises genomes into separate databases, so it is necessary
    to obtain the database from the associated metadata for each genome, along with the NCBI
    taxonomy id, which is need by kraken.

    The genome name is retrieved from the 'name' field of the 'org' element of the genome record.
    The taxonomy id is retrieved from the 'taxonomy_id' field of the 'org' element of the genome record.
    The database is retrieved from the 'dbname' attribute of the 'databases' elemennt of the genome record.
    
    required parameters:
        metadata_file: name of metadata file to parse

    returns:
        genome_dat(dict): dict keyed on genome name, which contains anonymous dicts 
        with 'name' (dbname) and 'taxid' keys
    """
    print('Parsing metadata....')
    infile=open('{}/{}'.format(data_dir,metadata_file))
    metadata=infile.read()
    doc = BeautifulSoup(metadata,'xml')
    genomes=doc.find_all('genome')

    genome_dat={}

    for genome in genomes:
        db=genome.find_all('databases')
        org=genome.find('organism')
        taxid=org['taxonomy_id']

        for d in db:
            if d['type']=='core':
                dbname=d['dbname']
                dbname = re.sub(r'_core_.*','',dbname, re.S)

            genome_dat[org['name']]={'name': dbname, 'taxid': taxid}

    pickle_file=metadata_file.replace('xml','pkl')
    pickle.dump( genome_dat, open( "{}/{}".format(data_dir,pickle_file), "wb" ) )

    return(genome_dat)
        
    
def download_file(ftphost, ftpdir, ftpfile):
    
    """
    Download file via ftp if newer than existing file, storing remote timestamp
    in filename.timestamp file
    
    arguments:
        ftphost: name of ftp server 
        ftpdir: path to directory of interest on ftp server
        ftpfile: Name of file to retreive
    
    returns:
    none
    """
    
    ftp = FTP(ftphost)
    ftp.login()
    ftp.cwd(ftpdir)
    timestamp = ftp.sendcmd('MDTM '+ ftpfile)
    # Successful command return will be prefixed with '213' return code
    if timestamp.startswith('213 '):
        timestamp=timestamp.replace('213 ', '')
    else:
        print("Error retrieving timestamp of %s" % ftpfile)
        exit(1)
    
    ts_file = ('%s/%s.timestamp' % (data_dir,ftpfile))
    if os.path.isfile(ts_file):
        with open(ts_file, 'r') as ts:
            old_ts = ts.read()
    else:
        old_ts = ''
            
    if timestamp != old_ts:
        print("Downloading new version of %s..." % ftpfile )
        outfile = open('%s/%s' % (data_dir,ftpfile), 'wb')
            
        ftp.retrbinary('RETR %s' % ftpfile, outfile.write)
        outfile.close()
            
        # save timestamp of downloaded file for future comparions
        ts=open(ts_file,'w')
        ts.write(timestamp)
        ts.close()
        
    
def download_metadata(metadata_file,div):
    
    """
    Download ensembl metadata (if newer than local version)
    
    arguments:
    metadata_file - name of metadata xml file
    div -- division to query
    
    returns:
    None
    """
    
    download_file('ftp.ensemblgenomes.org', '/pub/%s/current' % div, metadata_file)
        
    return()

def get_ensembl_release(div):
    
    """
    Retrieve ensembl release number from 'current' symlink on ftp server
    
    arguments:
    div -- division to query
    
    returns:
    release - string
    """
    
    # The ensembl release version can be found from the path on the FTP served
    # symlinked to the 'current' directory...
    
    ftp = FTP('ftp.ensemblgenomes.org')
    ftp.login()
    ftp.cwd('/pub/%s/current' % div)
    path=(ftp.pwd())
    
    regexp = r"release-([\d]+)"
    match =  re.search(regexp, path)
    release = match.group(1)
    
    return(release)

def add_taxid_to_records(genome,taxid):

    '''
    Update fasta headers of sequence records to incorporate taxid in kraken-friendly manner
    '''

    genome_file=os.path.basename(genome)
    genome_file=genome_file.replace('.gz','')
    sect=genome_file[:1]

    if not os.path.exists('{}/ensembl/{}'.format(data_dir,sect)):
        os.mkdir('{}/ensembl/{}'.format(data_dir,sect))
    
    out_handle=open('{}/ensembl/{}/{}'.format(data_dir,sect,genome_file),'w')

    with gzip.open(genome,'rt') as in_handle:
        for record in SeqIO.parse(in_handle, 'fasta'):
            record.id='{}|kraken:taxid|{}'.format(record.id,taxid)
            SeqIO.write(record,out_handle,'fasta')

    os.remove(genome)

def download_genome(genome, db, taxid, div):
    
    '''
    Download specified genome from ensembl...
    
    none
    '''

    # flag is created following successful download to allow us to resume 
    # failed attempts...
    flag='{}/flags/{}'.format(data_dir,genome)
    if not os.path.exists(flag):
        ftp_dir='/pub/%s/current/fasta/' % (div)
        ftp = FTP('ftp.ensemblgenomes.org')
        for x in range(retryCount):
            try:
                if x>0:
                    print('Login retry: attempt {}'.format(x))
                ftp.login()
                break
            except ftplib.all_errors:
                time.sleep(retrySleep)
        
        pat = r'.*.dna.toplevel.fa.gz$'
        
        if db==genome:
            species_dir = ('%s/%s/dna' % (ftp_dir,db))
        else:
            species_dir = ('%s/%s/%s/dna' % (ftp_dir,db,genome))
        
        ftp_files=None	
        for x in range(retryCount):
            try:
                if x>0:
                    print('nlst retry: attempt {}'.format(x))
                ftp_files=ftp.nlst(species_dir)
            except ftplib.all_errors:
                time.sleep(retrySleep)
            else:
                break

        dl_file=''
        if ftp_files is not None:
            for f in ftp_files:
                if re.match(pat,f):
                    dl_file = f
                    break

        if dl_file == '':
            print('dl_file not found (genome=%s; db=%s)' % (genome,db))
            pprint.pprint(ftp_files)
            
        ftp_file=ntpath.basename(dl_file)
        if ftp_file:    
            sect=genome[:1] 
            if not os.path.exists('%s/tmp/%s' % (data_dir,sect)):
                os.makedirs('%s/tmp/%s' % (data_dir,sect))

            ensembl_file=ftp_file.replace('.gz','')
            if not os.path.exists('%s/ensembl/%s/%s' % (data_dir,sect,ensembl_file)):
                fn='%s/tmp/%s/%s' % (data_dir,sect,ftp_file)
                outfile = open(fn, 'wb')
                for x in range(retryCount):
                    try:
                        if x>0:
                            print('Download retry ({}): attempt {}'.format(dl_file,x))
                        ftp.retrbinary('RETR %s' % dl_file, outfile.write)
                        break
                    except ftplib.all_errors as e:
                        print('Download error' + str(e))
                        time.sleep(retrySleep)
            
                outfile.close()
                add_taxid_to_records(fn,taxid)

            else:
                print('No ftp file(genome={}, db={})'.format(genome,db))
    Path(flag).touch()

def cleanup():

    shutil.rmtree('{}/ensembl'.format(data_dir))
    shutil.rmtree('{}/tmp'.format(data_dir))
    sys.exit(0)

if __name__ == "__main__":
    main()
    
