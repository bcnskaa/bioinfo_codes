#########################
#
# Copyright (c) 2013, SK Woolf <bcnskaa@gmail.com>.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#########################
import sys
import _thread
import time
import traceback
import math
import glob
from os import listdir
import os



# Summary:
#   Parse SQL schema and update UCSC BED data with the schema
#
# Input:
#   1) SQL schema or .TBL file from UCSC genome browser 
#   2) Corresponding BED file for the schema
# Output:
#   An annotated BED file
#
# Usage: 
#   Call this script 
#   >import sql_to_table
#   >sql_to_table.append_header_to_bed("cosmicRaw.sql", "cosmicRaw.txt", "cosmicRaw.txt2")

def main(argv):
    #do_scan(argv[0])
    #print("Arg: ", len(argv), argv)
    
    if(len(argv) == 3):
        append_header_to_bed(argv[0],argv[1],argv[2])
    elif len(argv) == 1 and ('.sql' in argv[0] or '.tbl' in argv[0]):
        #file_pattern = argv[0].replace("*", "")
        file_pattern = argv[0]
        # Scan through current directory
        schema_files = [f for f in os.listdir('.') if file_pattern in f]
        
        schema_names = [name.replace(file_pattern, "") for name in schema_files]
        
        for schema_name in schema_names:
            append_header_to_bed(schema_name + file_pattern, schema_name + ".txt", schema_name + ".txt2")
    else:
        print("Usage: SQL(.sql)|TABLE(.tbl) TAB-TXT OUT-TXT")


def read_sql_schema(schema_sql_fn):
    print("Reading from SQL schema ", schema_sql_fn, "...")
    
    _sql_schema = open(schema_sql_fn).read()
    _sql_schema = _sql_schema.splitlines()
    
    print("Number of line: ", len(_sql_schema))

    state = 0
    columns = []
    
    # Looking for specific keywords: KEY, CREATE TABLE, ` 
    for i in range(0, len(_sql_schema)):
        # Reach the end 
        if _sql_schema[i].startswith('  KEY '):
            state = 2
        
        # Other scenerio for breaking the loop 
        if state == 1 and not(_sql_schema[i].startswith('  `')):
            state = 2    
        
        # Break the loop
        if state == 2:
            break
        
        # Grab the name
        if state == 1:
            spos = len("  `")
            epos = _sql_schema[i].index("`", spos + 1)
            columns.append(_sql_schema[i][spos:epos])
        
        # Start of a table schema
        if _sql_schema[i].startswith('CREATE TABLE'):
            state = 1

            
    return columns


#
# tbl format can be prepared by copying the table schema from UCSC table browser
#
def read_tbl_schema(schema_tbl_fn):
    print("Reading from TBL schema ", schema_tbl_fn, "...")
    
    _tbl_schema = open(schema_tbl_fn).read()
    _tbl_schema = _tbl_schema.splitlines()
    
    columns = []
    for i in range(0, len(_tbl_schema)):
        if i > 0:
            fields = _tbl_schema[i].split('\t')
            columns.append(fields[0])
    
    return columns
    
# 
def append_header_to_bed(header_schema_fn, data_in_fn, data_out_fn):
    if '.sql' in header_schema_fn:
        headers = read_sql_schema(header_schema_fn)
    else:
        headers = read_tbl_schema(header_schema_fn)  
    
    if len(headers) == 0:
        return
    
    print('Number of header column: ', len(headers))

    headers = "\t".join(headers)
    line_n = 0

    with open(data_in_fn, "r") as f_in, open(data_out_fn, "w") as f_out:
        print('Reading data', data_in_fn, "and export to", data_out_fn)
        f_out.write(headers + "\n")    
        
        for line in f_in:
            f_out.write(line)
            line_n += 1


# Scan the directory dir and look for files ended with .sql
def do_scan(dir):
    names = [f.replace('.sql', '') for f in listdir(dir) if f.endswith('.sql')]
    for i in range(len(names)):
        print('processing ', n, '...')
        n = names[i]
        sql_fn = n + '.sql'
        data_fn = n + '.txt'
        out_fn = n + '.txt2'
        append_header_to_bed(sql_fn, data_fn, out_fn)


#import sql_to_table
if __name__ == "__main__":
    main(sys.argv[1:])
