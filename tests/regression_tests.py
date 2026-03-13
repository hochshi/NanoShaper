from ns_test import ns_test
import sys 
import argparse
import os
import numpy as np


def argsParser():
    parser = argparse.ArgumentParser(prog='nanoshaper_regression_tests', description='Manage regression tests for nanoshaper')
    parser.add_argument('--testname', help='specific regression test to run', default='all')
    parser.add_argument('--refdir', help='path of the reference output files/directories', default='./')
    parser.add_argument('--testdir', help='path of the output files/directories of the tested exec', default='./')
    parser.add_argument('--execdir', help='path of the executable file', default='../build/')
    parser.add_argument('--testtype', help='light test or full test', default='light')
    parser.add_argument('--remove', help='remove files outputted by the exec to test, or not', default='yes')
    parser.add_argument('--verbose', help='output timings and memory spaces, or not', default='no')
    args = parser.parse_args()

    return args


def runTest(test, ref_dir, test_dir, exec_dir, remove_output_test_files, verbose):
    
    """
    init a new ns_test obj
    run the given test

    :param test name of the test to be performed
    :param ref_dir absolute path of the directory where the reference tests are
    :param test_dir absolute path of the directory where the tests are to be compared with the reference tests
    :param exec_dir absolute path wherein the NanoShaper executable is available
    """

    _confprm = test + '_conf.prm'
    _ns_test = ns_test(_confprm, test, ref_dir, test_dir, exec_dir, remove_output_test_files, verbose)
    _ns_test.run()

    return _ns_test


if __name__ == "__main__":
    
    _ns_test_list = []

    args = argsParser()
    ref_dir = os.path.abspath(args.refdir)
    mother_test_dir = os.path.abspath(args.testdir)
    exec_dir = os.path.abspath(args.execdir)
    test_type = args.testtype

    if (args.remove == 'yes'):
        remove_output_test_files = True
    else:
        remove_output_test_files = False

    if (args.verbose == 'yes'):
        verbose = True
    else:
        verbose = False
    
    print("\nNanoShaper Regression Tests \n")

    try:
        if args.testname == 'all':
            for test in os.listdir(ref_dir):
                if test == '__pycache__':
                    continue
                test_dir = os.path.join(mother_test_dir, test)
                local_ref_dir = os.path.join(ref_dir, test)

                if os.path.isdir(test_dir):
                    _ns_test = runTest(test, local_ref_dir, test_dir, exec_dir, remove_output_test_files, verbose)
                    _ns_test_list.append(_ns_test)

        else:
            test = args.testname
            test_dir = os.path.join(mother_test_dir, test)
            ref_dir = os.path.join(ref_dir, test)
            print("ref and test dirs: ")
            print(ref_dir, test_dir)

            if not os.path.isdir(test_dir):
                raise Exception("The test directory does not exist\n") 
            else:
                _ns_test = runTest(test, ref_dir, test_dir, exec_dir, remove_output_test_files, verbose)
                _ns_test_list.append(_ns_test)

        os.chdir(mother_test_dir)
        N = len(_ns_test_list)
        failed_counter = 0

        for i_ns_test in _ns_test_list:
            if i_ns_test.failed:
                failed_counter = failed_counter + 1

        print("Regression tests passed: {} / {}\n".format((N - failed_counter), N))

    except Exception as err:
        print("regression_tests.py failed with the following error: {}".format(err))
        sys.exit(-1)

   
    if failed_counter > 0:
        raise Exception("Some regression tests failed!!!!") 



    
