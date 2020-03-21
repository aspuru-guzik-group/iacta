import numpy as np
import subprocess
import shutil
import os
import tempfile

class xtb_run:
    def __init__(self,xtb, geom_file, *args, cwd=None, xcontrol=None):
        self.dir = tempfile.mkdtemp(dir=cwd)
        self.out = open(self.dir + "/xtb.out", "w")
        self.err = open(self.dir + "/xtb.err", "w")
        self.coord = shutil.copy(geom_file, self.dir)
        
        self.args = [xtb, os.path.basename(self.coord)]
        if xcontrol:
            self.xcontrol = shutil.copy(xcontrol, self.dir)
            self.args += ["-I", os.path.basename(self.xcontrol)]
        else:
            self.xcontrol = None

        self.args += args
        
        
        self.proc = subprocess.Popen(
            self.args,
            stderr=self.err,
            stdout=self.out,
            cwd=self.dir)nn

    def close(self):
        self.out.close()
        self.err.close()
        shutil.rmtree(self.dir)


class xtb_driver:
    def __init__(self, path_to_xtb):
        self.defaults = ["-T", "1"]
        self.xtb_bin = path_to_xtb

    def add_to_defaults(args):
        self.defaults.append(args)


