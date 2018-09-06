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
@click.option('--skip-phase',default=False,is_flag=True,help='Flag to skip phasing')
def create(name,vcf,fasta,skip_phase):
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
        h = hap.HapDab.from_files(name,vcf,fasta,skip_phase=skip_phase)
cli.add_command(create)

#----------------------------
#    Impute Commands
#----------------------------
@click.command(short_help='Impute a VCF file.')
@click.argument('name',metavar='<name>')
@click.argument('vcf',metavar='<vcf>')
def impute(name,vcf):
    '''
    Impute a VCF file to a higher genotype density.

    \b
    Positional Arguments:
    <name> - The name of the HapDab databse to use for imputation.
    <vcf>  - The input (low density) VCF to be imputed
    '''
    pass

cli.add_command(impute)

#----------------------------
#    Remap Commands
#----------------------------
@click.command(help='Remap SNP coordinates.')
#@click.argument('dtype',metavar='<dtype>')
#@click.argument('name',metavar='<name>')
def remap(name, dtype):
    pass

cli.add_command(remap)


