import numpy as np
import subprocess
import shutil
import os
import tempfile


def make_xcontrol(xcontrol_dictionary, fn):
    f = open(fn, "w")
    for key, val in xcontrol_dictionary.items():
        if type(val) == tuple:
            f.write("$" + str(key) + "\n")
            for v in val:
                f.write(v + "\n")
            f.write("$end\n")
        else:
            f.write("$" + str(key) + str(val) + "\n")

    f.close()
    return fn
            

class xtb_run:
    def __init__(self, xtb, geom_file,
                 *args, scratch=".", cwd=".", xcontrol=None,
                 return_files=[], before_geometry="--",
                 block=True):
        self.dir = tempfile.mkdtemp(dir=scratch)
        self.out = open(self.dir + "/xtb.out", "w")
        self.err = open(self.dir + "/xtb.err", "w")
        self.coord = shutil.copy(geom_file, self.dir)
        self.cwd = cwd
        self.return_files = return_files   # list of files to take out when
                                           # run finishes
        
        self.args = [xtb]
        if xcontrol:
            self.xcontrol = make_xcontrol(xcontrol, self.dir + "/.xcontrol")
            self.args += ["-I", os.path.basename(self.xcontrol)]
        else:
            self.xcontrol = None

        self.args += args
        self.args += [before_geometry, os.path.basename(self.coord)]
        
        self.kwargs = dict(stderr=self.err,
                           stdout=self.out,
                           cwd=self.dir)
        
        self.proc = None

    def start(self, blocking=True):
        assert self.proc is None
        self.proc = subprocess.Popen(self.args, **self.kwargs)
        if blocking:
            self.proc.wait()

    def kill(self):
        self.assert_running()
        self.proc.kill()

    def assert_done(self):
        if self.proc is None:
            raise RuntimeError("Process not started")
        
        if self.proc.poll() is None:
            raise RuntimeError("Process not finished")

    def assert_running(self):
        if self.proc is None:
            raise RuntimeError("Process not started")
        
        if self.proc.poll() is None:
            pass
        else:
            raise RuntimeError("Process finished")

    def cp(self, file_in, file_out=""):
        # Copy a file from xtb to file_out in the parent directory.
        self.assert_done()
        return shutil.copy(self.dir + "/" + file_in,
                           self.cwd + "/" + file_out)

    def close(self):
        self.assert_done()
        output = []
        for file_in, file_out in self.return_files:
            output += [self.cp(file_in,file_out)]
        
        self.out.close()
        self.err.close()
        shutil.rmtree(self.dir)
        return output

    def __call__(self, blocking=True):
        # Run xtb in a blocking manner, extract return files, close temporary
        # directories.
        self.start(blocking=blocking)
        output = self.close()
        
        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)



class xtb_driver:
    def __init__(self, path_to_xtb="", xtb_args=[],
                 scratch="."):
        self.extra_args = xtb_args
        self.xtb_bin = path_to_xtb + "xtb"
        self.crest_bin = path_to_xtb + "crest"
        self.scratchdir = scratch

    def optimize(self, geom_file, out_file,
                 xcontrol=None,
                 level=None,
                 compute_hessian=False,
                 log=None):
        
        file_ext = geom_file[-3:]
        return_files=[("xtbopt." + file_ext, out_file)]
        if log:
            return_files += [("xtbopt.log", log)]

        if compute_hessian:
            oflag = "--ohess"
            if level is None:
                level="vtight"

            return_files += [("hessian", "hessian_" + out_file)]
        else:
            oflag = "--opt"
            if level is None:
                level="normal"
            
        opt = xtb_run(self.xtb_bin, geom_file, 
                      oflag, level,
                      *self.extra_args,
                      xcontrol=xcontrol,
                      scratch=self.scratchdir,
                      return_files=return_files)
        return opt

    def multi_optimize(self, geom_file, out_file,
                       xcontrol=None,
                       level=None):
        file_ext = geom_file[-3:]
        return_files=[("crest_ensemble." + file_ext, out_file)]

        if level is None:
            level="normal"
            
        opt = xtb_run(self.crest_bin,
                      geom_file,
                      "-opt", level,
                      *self.extra_args,
                      before_geometry="-mdopt",
                      xcontrol=xcontrol,
                      scratch=self.scratchdir,                      
                      return_files=return_files)
        return opt

    def metadyn(self, geom_file, out_file, xcontrol=None):
        return_files=[("xtb.trj", out_file)]
        md = xtb_run("xtb", geom_file,
                     "--metadyn",
                     xcontrol=xcontrol,
                     scratch=self.scratchdir,                     
                     return_files=return_files)
        return md
