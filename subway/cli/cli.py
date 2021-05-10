import click
from rich import print


@click.command()
def entry_point():
    print(":metro: Subway CLI!")
