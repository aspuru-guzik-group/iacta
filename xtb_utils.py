import numpy as np
import subprocess
import shutil
import os
import tempfile

class xtb_run:
    def __init__(self,xtb, geom_file, *args, cwd=".", xcontrol=None):
        self.dir = tempfile.mkdtemp(dir=cwd)
        self.out = open(self.dir + "/xtb.out", "w")
        self.err = open(self.dir + "/xtb.err", "w")
        self.coord = shutil.copy(geom_file, self.dir)
        self.cwd = cwd
        
        self.args = [xtb]
        if xcontrol:
            self.xcontrol = shutil.copy(xcontrol, self.dir)
            self.args += ["-I", os.path.basename(self.xcontrol)]
        else:
            self.xcontrol = None

        self.args += args
        self.args += ["--", os.path.basename(self.coord)]
        
        
        self.proc = subprocess.Popen(
            self.args,
            stderr=self.err,
            stdout=self.out,
            cwd=self.dir)

    def cp(self, file_in, file_out=""):
        # Copy a file from xtb to file_out in the parent directory.
        return shutil.copy(self.dir + "/" + file_in,
                    self.cwd + "/" + file_out)


    def close(self):
        self.out.close()
        self.err.close()
        shutil.rmtree(self.dir)


class xtb_driver:
    def __init__(self, path_to_xtb="xtb", xtb_args=["-T 1"]):
        self.extra_args = xtb_args
        self.xtb_bin = path_to_xtb

    def optimize(self, geom_file, out_file,
                 xcontrol=None,
                 log=None):
        file_ext = geom_file[-3:]
        opt = xtb_run(self.xtb_bin, geom_file,
                      "--opt", xcontrol=xcontrol,
                      *self.extra_args)
        opt.proc.wait()
        out = opt.cp("xtbopt." + file_ext, out_file)
        if log:
            opt.cp("xtbopt.log", log)
        opt.close()
        return out 


