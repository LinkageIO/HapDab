#!/usr/bin/env python3

import click
import shutil
import hapdab as hap
import minus80 as m80


@click.group(epilog=f'Version {hap.__version__}\n{hap.__file__}')
def cli():
    '''
    \b
         _                     _       _     
        | |__   __ _ _ __   __| | __ _| |__  
        | '_ \ / _` | '_ \ / _` |/ _` | '_ \ 
        | | | | (_| | |_) | (_| | (_| | |_) |
        |_| |_|\__,_| .__/ \__,_|\__,_|_.__/ 
                    |_|   (powered by minus80) 
    '''


#----------------------------
#    create Commands
#----------------------------
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
    try:
        if name is None or vcf is None or fasta is None:
            click.echo('Invalid Syntax: use --help for more.')
        else:
            h = hap.HapDab.from_files(name,vcf,fasta,skip_phase=skip_phase)
    except Exception as e:
        import ipdb; ipdb.set_trace()
cli.add_command(create)

#----------------------------
#    List Commands
#----------------------------
@click.command(short_help='List the available hapdab datasets',
    help='Reports the available **Frozen** datasets.'
)
@click.option('--name',  default=None,
    help="The name of the dataset you want to check is available. The default value is the wildcard '*' which will return all available datasets with the specified dtype."
)
def list(name):
    m80.Tools.available(dtype='HapDab',name=name)
cli.add_command(list)



#----------------------------
#    Impute Commands
#----------------------------
@click.command(short_help='Impute a VCF file.')
@click.argument('name',metavar='<name>')
@click.argument('vcf',metavar='<vcf>')
@click.option('--out', default=None,
    help=('An optional output name. If not specified, one '
           'will be created from <vcf>')
)
def impute(name,vcf,out):
    '''
    Impute a VCF file to a higher genotype density.

    \b
    Positional Arguments:
    <name> - The name of the HapDab databse to use for imputation.
    <vcf>  - The input (low density) VCF to be imputed
    '''
    h = hap.HapDab(name)
    if out is None:
        out = vcf.replace('.vcf','').replace('.gz','') + '.imputed.vcf.gz'
    imputed = h.impute(vcf)
    shutil.copyfile(imputed,out)


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


#----------------------------
#    Cloud Commands
#----------------------------
@click.group()
def cloud():
    '''
    Manage your HapDab datasets in the cloud.
    '''

cli.add_command(cloud)

@click.command()
@click.option('--name', metavar='<name>',default=None)
def list(name):
    '''List available datasets'''
    cloud = m80.CloudData()
    cloud.list(
        dtype='HapDab',
        name=name,
        raw=False
    )


@click.command()
@click.argument('name', metavar='<name>')
def push(name):
    '''
    \b
    Push a HapDab dataset to the cloud.

    \b
    Positional Arguments:
    <name> - the name of the m80 dataset or raw filename if --raw is set.
    '''
    cloud = m80.CloudData()
    cloud.push(
        'HapDab',
        name,
        raw=False
    )

@click.command()
@click.argument('name', metavar='<name>')
@click.option('--output',default=None,help="Output filename, defaults to <name>. Only valid with --raw")
def pull(name,output):
    '''
    Pull a HapDab dataset from the cloud.
    '''
    cloud = m80.CloudData()
    cloud.pull(
        'HapDab',
        name,
        raw=False,
        output=output
    )

@click.command()
@click.argument('name', metavar='<name>')
def remove(name):
    '''
    Delete a HapDab dataset from the cloud.
    '''
    cloud = m80.CloudData()
    cloud.remove(
        'HapDab',
        name,
        raw=False
    )


cloud.add_command(list)
cloud.add_command(push)
cloud.add_command(pull)
cloud.add_command(remove)

cli.add_command(cloud)

