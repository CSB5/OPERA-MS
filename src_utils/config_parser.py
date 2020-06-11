import argparse

parser = argparse.ArgumentParser(add_help=False)

mandatory = parser.add_argument_group("mandatory arguments")
mandatory.add_argument("-c", "--config",
                       #'config',
                       required=True,
                       help="config file")
