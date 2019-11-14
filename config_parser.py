import argparse

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-t", "--thread", help="Number of thread [default value in config file]")
#parser.add_argument("-d", "--thread", help="Number of thread [default value in config file]")
mandatory = parser.add_argument_group("mandatory arguments")
mandatory.add_argument("-c", "--config",
                       #'config',
                       required=True,
                       help="config file")
