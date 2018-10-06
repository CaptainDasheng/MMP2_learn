#!/usr/bin/env python  

from optparse import OptionParser, OptionGroup

def testfuc1():

    
    parser = OptionParser("usage: %jump2 -r [options] args")
    group1 = OptionGroup(parser, 'required flags')

    parser_ext = OptionParser("usage: %jump2 -c [options] args")
    group2 = OptionGroup(parser, 'required flags')
    parser_chk = OptionParser("usage: %jump2 -k [options] args")
    group3 = OptionGroup(parser, 'required flags')

    group2.add_option('-j', '--job', dest='jobs', action='store_true',default=False, help='test help')
    group3.add_option('-g', '--gap', dest='gap', action='store_true', default=False, help='test help')
    group1.add_option('-l', '--launch', dest='sumbit', action='store_true', default=False, help='test help')

    parser.add_option_group(group1)
    parser.add_option_group(group2)
    parser.add_option_group(group3)
    (options,args) = parser.parse_args()

testfuc1()


