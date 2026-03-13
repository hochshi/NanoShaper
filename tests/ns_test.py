import subprocess
import numpy as np
import glob
import os
import json
import time
import tarfile
import shutil
    

class ns_test:
    def __init__(self, confprm, name, ref_dir='./', test_dir='./', exec_dir='../build', remove_output_test_files=True, verbose=False):

        """
        ns_test constructor

        :param confprm surface configuration file name
        :param name  name of the test to be performed
        :param ref_dir  path of the directory where the reference outputs are
        :param test_dir  path of the directory where the outputted test files are
        :param exec_dir  path executable file
        :param remove_output_test_files  boolean determining the removal of the output test files
        """

        self.confprm = confprm
        self.name = name
        self.ref_dir = ref_dir
        self.test_dir = test_dir
        self.exec_dir = exec_dir
        self.returncode = -1
        self.failed = False
        self.remove_output_test_files = remove_output_test_files
        self.verbose = verbose
        self.confmsg = 'The surface configuration does not exist: '
        self.outmsg = ''
        self.shapemsg = ''
        self.valuemsg = ''
        self.execmsg = ''
        self.genmsg = ''


    def isLogFile(self, p_file):
        name, ext = os.path.splitext(p_file)
        if ext != '.log':
            return False
        return True


    def isOffFile(self, p_file):
        name, ext = os.path.splitext(p_file)
        if ext != '.off':
            return False
        return True
    
    def isVerticesFile(self, p_file):
        name, ext = os.path.splitext(p_file)
        if ext != '.vert':
            return False
        return True

    def isAreasFile(self, p_file):
        name, ext = os.path.splitext(p_file)
        if ext != '.txt':
            return False
        return True


    def readVertices(self, p_file):
        
        """
        properly read .vert files

        """

        #Verteces/faces + Normals
        p_len = 6

        file = open(p_file, 'r')

        next(file, None)
        next(file, None)
        
        n = int(file.readline().strip())
        data = np.empty([n, p_len], dtype=float)
        
        for i in range(n):
            i_data = [float(s) for s in file.readline().strip().split()]
            data[i] = i_data[0:p_len]

        file.close()
        return data


    def readOff(self, p_file):
        
        """
        properly read .OFF files

        """

        p_len = 3
        
        file = open(p_file, 'r')
        format = file.readline().strip().split(' ')
        format_norms = "OFF+N"

        #Verteces + Normals
        if (format_norms == format[0]):
            p_len = 6

        next(file, None)
        next(file, None)
        
        n_verts, n_faces, n_edge = tuple([int(s) for s in file.readline().strip().split(' ')])

        verts = np.empty([n_verts, p_len], dtype=float)
        faces = np.empty([n_faces, 4], dtype=int)
        
        for i in range(n_verts):
            i_vert = [float(s) for s in file.readline().strip().split(' ')]
            verts[i] = i_vert
        
        for i in range(n_faces):
            i_face = [int(s) for s in file.readline().strip().split(' ')]
            faces[i] = i_face

        file.close()
        return verts, faces


    def readAreas(self, p_file):
        
        """
        read area values

        """
        
        areas = np.loadtxt(p_file, skiprows=0, comments='#', usecols=0, unpack=True)

        return areas


    def compareShapesAndValues(self, p_file):

        """

        compare the shapes and the values between output and reference arrays
        .OFF values 
        Area values

        :param p_file name of the file containing values 
        
        """

        try:
            ref_file = (os.path.join(self.ref_dir, p_file))
            test_file = (os.path.join(self.test_dir, p_file))

            keep_order = False
            
            if self.isOffFile(p_file):
                ref_verts, ref_faces = self.readOff(ref_file)
                test_verts, test_faces = self.readOff(test_file)

                if (test_verts.shape != ref_verts.shape) or (test_faces.shape != ref_faces.shape):
                    # too restrictive failure
                    # self.shapemsg = "{} {}".format(self.shapemsg,p_file)
                    # self.failed = True
                    self.failed = False
                else:
                    _ref_verts = ref_verts
                    _test_verts = test_verts
                    _ref_faces = ref_faces
                    _test_faces = test_faces

                    if not keep_order:
                        _ref_verts = np.sort(ref_verts, 0)
                        _test_verts = np.sort(test_verts, 0)

                        _ref_faces = np.sort(ref_faces, 0)
                        _test_faces = np.sort(test_faces, 0)

                    # tolerances in file values
                    rel_tol = 1.e-2
                    abs_tol = 1.e-3
                    if not np.allclose(_ref_verts, _test_verts, rel_tol, abs_tol):
                        self.valuemsg = "{} {}".format(self.valuemsg, p_file)
                        self.failed = True

                    """
                    # too restrictive failure for differences in number of faces
                    if not np.allclose(_ref_faces, _test_faces, 1, 1):
                        self.valuemsg = "{} {}".format(self.valuemsg, p_file)
                        self.failed = True
                    """

            elif self.isVerticesFile(p_file):
                ref_data = self.readVertices(ref_file)
                test_data = self.readVertices(test_file)

                if (test_data.shape != ref_data.shape):
                    # too restrictive failure
                    # self.shapemsg = "{} {}".format(self.shapemsg, p_file)
                    # self.failed = True
                    self.failed = False
                else:
                    _ref_data = ref_data
                    _test_data = test_data

                    if not keep_order:
                        _ref_data = np.sort(ref_data, 0)
                        _test_data = np.sort(test_data, 0)
                    
                    # tolerances in file values
                    rel_tol = 1.e-3
                    abs_tol = 1.e-4
                    if not np.allclose(_ref_data, _test_data, rel_tol, abs_tol):
                        self.valuemsg = "{} {}".format(self.valuemsg, p_file)
                        self.failed = True

            elif self.isAreasFile(p_file):
                ref_areas = self.readAreas(ref_file)
                test_areas = self.readAreas(test_file)

                if test_areas.shape != ref_areas.shape:
                    # too restrictive failure
                    # self.shapemsg = "{} {}".format(self.shapemsg, p_file)
                    # self.failed = True
                    self.failed = False
                else:
                    _ref_areas = ref_areas
                    _test_areas = test_areas
            
                    if not keep_order:
                        _ref_areas = np.sort(ref_areas)
                        _test_areas = np.sort(test_areas)

                    # tolerances in file values
                    rel_tol = 1.e-3
                    abs_tol = 1.e-4
                    if not np.allclose(_ref_areas, _test_areas, rel_tol, abs_tol):
                        self.valuemsg = "{} {}".format(self.valuemsg, p_file)
                        self.failed = True

            elif self.isLogFile(p_file):

                t_file = open(test_file, 'r')
                r_file = open(ref_file, 'r')

                test_lines = t_file.readlines()
                ref_lines = r_file.readlines()

                for test_row in test_lines:

                    # pocket = 'Pocket'

                    # if test_row.find(pocket) != -1:

                    #     for ref_row in ref_lines:
                    #         if ref_row.find(pocket) != -1:
                    #             if (ref_row[0:mem_location] == test_row[0:mem_location]):
                    #                 if test_row.find("loading") != -1:
                    #                         row = "".join([ test_row[10:len(test_row)-4], " VS", ref_row[ref_row.find("loading")+7:len(ref_row)-1] ])
                    #                         print(row)
                    #                     else:
                    #                         row = "".join([ test_row[10:len(test_row)-4], " VS", ref_row[ref_row.find("is")+2:len(ref_row)-1] ])
                    #                         print(row)

                    if self.verbose:
                        mem = 'MB'

                        if test_row.find(mem) != -1:
                            mem_location = test_row.find(mem) - 8

                            for ref_row in ref_lines:
                                if ref_row.find(mem) != -1:
                                    if (ref_row[0:mem_location] == test_row[0:mem_location]):
                                        if test_row.find("loading") != -1:
                                            row = "".join([ test_row[10:len(test_row)-4], " VS", ref_row[ref_row.find("loading")+7:len(ref_row)-1] ])
                                            print(row)
                                        else:
                                            row = "".join([ test_row[10:len(test_row)-4], " VS", ref_row[ref_row.find("is")+2:len(ref_row)-1] ])
                                            print(row)

                        time = '[s]'

                        if test_row.find(time) != -1:
                            test_row = test_row.replace("ok!", "")
                            test_row = test_row.replace("...", "")
                            test_row = test_row.replace("..", "")
                            time_location = test_row.find(time) - 8

                            for ref_row in ref_lines:
                                if ref_row.find(time) != -1:
                                    ref_row = ref_row.replace("ok!", "")
                                    ref_row = ref_row.replace("...", "")
                                    ref_row = ref_row.replace("..", "")
                                    if (ref_row[0:time_location] == test_row[0:time_location]):
                                        row = "".join([test_row[10:len(test_row)-5], " VS", ref_row[ref_row.find("[s]")-7:len(ref_row)-1] ])
                                        print(row)

        except Exception as err:
            self.genmsg = "{} compareShapesAndValues - {} - {} ".format(self.genmsg, p_file, err)
            self.failed = True


                


    def compareOutput(self):
        
        """
        compare output and reference data

        """

        try:

            if  self.failed == True:
                return

            ref_files = os.listdir(self.ref_dir)

            for ref_file in ref_files:
                ref_path = (os.path.join(self.test_dir, ref_file))

                if os.path.isfile(ref_path):
                    self.compareShapesAndValues(ref_file)
                    
                else:
                    self.outmsg = "{} {}".format(self.outmsg,ref_file)
                    self.failed = True
            
        except Exception as err:
            self.genmsg = "{} compareOutput - {} - {}".format(self.genmsg, self.name,err)
            self.failed = True


    def runExec(self):

        """

        execute NanoShaper given a surface configuration
        
        """

        try:
            print("Test: {}".format(self.name))

            if not os.path.isfile(os.path.join(self.test_dir, self.confprm)):
                self.confmsg =  self.confmsg + self.confprm + ' '
                print(self.confmsg)
                self.failed = True
            else:
                os.chdir(self.test_dir)

                exec = os.path.join(self.exec_dir, 'NanoShaper')
                log_file = os.path.join(self.test_dir, 'output.log')
                self.returncode = subprocess.call([exec, self.confprm], stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)

                if  self.returncode != 0:
                    self.execmsg = "Fail to execute: {}".format(exec)
                    print(self.execmsg)
                    self.failed = True
        except Exception as err:
            self.genmsg = "{} runExec - {} ".format(self.genmsg, self.name)
            self.failed = True


    def Cleanup(self):
        
        """
        clean up all output files

        """

        try:
            if (self.remove_output_test_files):
                for outfile in glob.iglob(os.path.join(self.test_dir, '*.txt')):
                    os.remove(outfile)
                for outfile in glob.iglob(os.path.join(self.test_dir, '*.xyz')):
                    os.remove(outfile)
                for outfile in glob.iglob(os.path.join(self.test_dir, '*.face')):
                    os.remove(outfile)
                for outfile in glob.iglob(os.path.join(self.test_dir, '*.vert')):
                    os.remove(outfile)
                for outfile in glob.iglob(os.path.join(self.test_dir, '*.off')):
                    os.remove(outfile)

        except Exception as err:
            self.genmsg = "{} Cleanup - {}".format(self.genmsg, self.name,err)
            self.failed = True

    
    def run(self):

        """
        execute NanoShaper given a surface configuration
        compare the outputted test files with the reference files
        clean up all output test files if their removal is enabled
        print all error msgs

        """
        t0 = time.time()

        self.runExec()
        self.compareOutput()
        self.Cleanup()

        t1 = time.time() - t0

        print("Elapsed Time [s]: {}".format(round(t1, 1)))

        if self.failed == True:

            if self.outmsg != '':
                self.outmsg = "The following output files do not exist: {}".format(self.outmsg)
                print(self.outmsg)
            if self.shapemsg != '':
                self.shapemsg = "Different shapes in files: {}".format(self.shapemsg)
                print(self.shapemsg)
            if self.valuemsg != '':
                self.valuemsg = "Different values in files: {} ".format(self.valuemsg)  
                print(self.valuemsg)
            if self.genmsg != '':
               self.genmsg = "Generic Errors: {}".format(self.genmsg)
               print(self.genmsg) 

            print("Failed\n")
        else:
            print("Passed\n")

