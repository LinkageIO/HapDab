#!/usr/bin/env python3

import click
import hapdab as hap


@click.group(epilog=f'Version {hap.__version__}\n{hap.__file__}')
def cli():
    '''
    \b
         _                     _       _     
        | |__   __ _ _ __   __| | __ _| |__  
        | '_ \ / _` | '_ \ / _` |/ _` | '_ \ 
        | | | | (_| | |_) | (_| | (_| | |_) |
        |_| |_|\__,_| .__/ \__,_|\__,_|_.__/ 
                    |_|    
    '''

@click.command(short_help='Create a HapDab database')
@click.argument('name',metavar='<name>')
@click.argument('vcf',metavar='<vcf>')
@click.argument('fasta',metavar='<fasta>')
def create(name,vcf,fasta):
    '''
    \b
    Create a HapDab database for imputation from
    a reference genome and a VCF file.

    \b
    Positional Arguments:
    <name>  - The name of the resultant HapDab database
    <vcf>   - The input VCF file
    <fasta> - An input Fasta file
    '''
    if name is None or vcf is None or fasta is None:
        click.echo('Invalid Syntax: use --help for more.')
    else:
        h = hap.HapDab.from_files(name,vcf,fasta)
cli.add_command(create)

#----------------------------
#    List Commands
#----------------------------
@click.command(short_help='Impute a VCF file.',
    help='Impute a VCF file to a higher genotype density.'
)
#@click.option('--name',  default=None,
#    help="The name of the dataset you want to check is available. The default value is the wildcard '*' which will return all available datasets with the specified dtype."
#)
#@click.option('--dtype', default=None,
#    help='Each dataset has a datatype associated with it. E.g.: `Cohort`. If no dtype is specified, all available dtypes  will be returned.'
#)
def impute():
    pass

cli.add_command(impute)

#----------------------------
#    delete Commands
#----------------------------
@click.command(help='Remap SNP coordinates.')
#@click.argument('dtype',metavar='<dtype>')
#@click.argument('name',metavar='<name>')
def remap(name, dtype):
    pass

cli.add_command(remap)


