import numpy as np
import subprocess
import shutil
import os
import tempfile

def make_xcontrol(xcontrol_dictionary, fn):
    """Transform a dictionary of parameters to an xTB xcontrol file.

    Parameters:
    -----------
    xcontrol_dictionary (dict) : dictionary of items to include in the
    xcontrol file. Values that are string will be rendered into $key val
    pairs, while values that are tuples will be rendered into blocks:
       $key
         val[1]
         val[2]
         ...
       $end

    fn (str) : filename to save xcontrol file to.

    Returns:
    --------
    str : filename fn.
    """

    
    f = open(fn, "w")
    for key, val in xcontrol_dictionary.items():
        if val is None:
            pass
        else:
            if type(val) == tuple:
                f.write("$" + str(key) + "\n")
                for v in val:
                    f.write(v + "\n")
                f.write("$end\n")
            else:
                f.write("$" + str(key) + " " + str(val) + "\n")

    f.close()
    return fn
            

class xtb_run:
    def __init__(self, xtb, geom_file, *args,
                 before_geometry="--",
                 scratch=".",
                 prefix="",
                 xcontrol=None,
                 restart=None,
                 failout=None,
                 delete=True,
                 return_files=[]):
        """Build a container for an xtb run.

        Note: the xtb job is only *prepared* when the object is defined. To run
        the job, use the __call__() or start(self) methods below.

        Parameters:
        ----------
        xtb (str) : Path to xtb or crest binary.

        geom_file (str) : Path to geometry file.

        *args : additional arguments to pass to xtb.

        Optional Parameters:
        --------------------
        before_geometry (str) : string placed before geom file, defaults to
        "--".

        scratch (str) : scratch directory where temporary files are placed.

        xcontrol (dict) : dictionary of xcontrol options, interpreted using
        make_xcontrol.

        failout (str) : Path where xtb.out is output if the run fails.

        delete (bool) : If true, delete all temporary files after run.
        Defaults to true

        return_files (list) : list of tuples of filenames (filein, fileout) of
        files to be copied out of the temporary run directory automatically
        when close() is called. For example, [("xtbopt.xyz", "my_opt.xyz")]
        will generate the file my_opt.xyz from the xtb optimized geometry
        xtbopt.xyz in the run directory.

        """

        self.dir = tempfile.mkdtemp(dir=scratch, prefix="tmp"+prefix)
        self.delete = delete
        self.coord = shutil.copy(geom_file, self.dir)
        self.failout = failout
        if restart:
            shutil.copy(restart, self.dir + "/xtbrestart")
            
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
        self.kwargs = dict(cwd=self.dir)

        # Holder for the process
        self.proc = None

    def start(self, blocking=True):
        """Start xtb job."""
        assert self.proc is None
        self.proc = subprocess.Popen(self.args, **self.kwargs)
        if blocking:
            self.proc.wait()

    def kill(self):
        """Kill currently running xtb job."""
        self.assert_running()
        self.proc.kill()

    def assert_done(self):
        """Raise RuntimeError if job is not done."""
        if self.proc is None:
            raise RuntimeError("Process not started")
        
        if self.proc.poll() is None:
            raise RuntimeError("Process not finished")

    def assert_running(self):
        """Raise RuntimeError if job is not currently running."""
        if self.proc is None:
            raise RuntimeError("Process not started")
        
        if self.proc.poll() is None:
            pass
        else:
            raise RuntimeError("Process finished")

    def cp(self, file_in, file_out=""):
        """Copy a file from the xtb temporary directory to file_out."""        
        self.assert_done()
        if file_out == "":
            file_out = os.path.basename(file_in)
        return shutil.copy(self.dir + "/" + file_in, file_out)

    def close(self, kill=False):
        """Cleanup xtb job."""
        if self.proc:        
            if kill:
                try:
                    self.kill()
                except RuntimeError:
                    pass
            else:
                self.assert_done()
            output = []
            for file_in, file_out in self.return_files:
                output += [self.cp(file_in,file_out)]

        if self.delete:
            shutil.rmtree(self.dir)
        else:
            self.xtb_log(self.dir + "/xtb.out",
                         self.dir + "/xtb.err")

        # Check error code != 0
        if self.proc.returncode and self.failout:
            self.xtb_log(self, self.failout)
            
        return output

    def xtb_log(self, fileout, fileerr=None):
        """Dump content of xtbout and xtberr to files."""
        # Check if xtb is finished
        self.assert_done()

        # Dump to file
        with open(fileout, "w") as f:
            f.write(self.proc.stdout)
            
        if fileerr:
            with open(fileerr, "w") as f:
                f.write(self.proc.stderr)
            

    def __call__(self, blocking=True):
        """Run xtb job and cleanup."""
        self.start(blocking=blocking)
        output = self.close()
        
        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)

class xtb_driver:
    def __init__(self, path_to_xtb_binaries="",
                 delete=True,
                 xtb_args=[], scratch="."):
        """Utility driver for xtb runs.

        Methods include various kind of xtb runs.

        Optional parameters:
        --------------------

        path_to_xtb (str) : absolute path to folder containing xtb and crest.
        Defaults to "" -> binaries are in $PATH.

        xtb_args (list of str): additional arguments to pass to xtb (for
        example, -gfn0, etc.). Default to [].

        delete (bool) : if true, delete run files after runs have completed.

        scratch (str) : scratch directory for xtb runs, defaults to ".".

        """
        self.extra_args = xtb_args
        self.xtb_bin = path_to_xtb_binaries + "xtb"
        self.crest_bin = path_to_xtb_binaries + "crest"
        self.scratchdir = scratch
        self.delete=delete

    def optimize(self,
                 geom_file,
                 out_file,
                 xcontrol=None,
                 level=None,
                 compute_hessian=False,
                 log=None,
                 failout=None,
                 restart=None):
        """Optimize a molecule.

        Parameters:
        -----------

        geom_file (str) : path to the file containing the molecular geometry.

        out_file (str): path to file where optimized geometry is saved.

        Optional Parameters:
        --------------------

        xcontrol (dict) : xcontrol dictionary to be interpreted by
        make_xcontrol.

        level (str) : Optimization level. Defaults to "normal" if
        compute_hessian is False, or to "tight" otherwise.

        compute_hessian (bool) : Compute Hessian if True. Hessian file is
        saved to hessian_+ out_file.

        log (str) : If set, additionally output the xtb geometry log.

        failout (str) : Path where xtb.out is output if the optimization
        fails.

        restart (str) : Setup a xtbrestart file that will be read for the run
        and written with restart details when the run is over.

        Returns:
        --------

        xtb_run : The optimization job. Run using xtb_run().

        """
        
        file_ext = geom_file[-3:]
        return_files=[("xtbopt." + file_ext, out_file)]
        if log:
            return_files += [("xtbopt.log", log)]
        if restart:
            return_files += [("xtbrestart", restart)]
            
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
                      restart=restart,
                      prefix="OPT",
                      scratch=self.scratchdir,
                      delete=self.delete,
                      failout=failout,
                      return_files=return_files)
        return opt

    def metadyn(self,
                geom_file, out_file,
                failout=None,
                xcontrol=None):
        """Run metadynamics on a molecule.

        Parameters:
        -----------

        geom_file (str) : path to the file containing the initial molecular
        geometries.

        out_file (str): path to file where results are saved.

        Optional Parameters:
        --------------------

        xcontrol (dict) : xcontrol dictionary to be interpreted by
        make_xcontrol.

        failout (str) : Path where xtb.out is output if the metadynamics
        fails.

        Returns:
        --------

        xtb_run : The metadynamics job. Run using xtb_run().

        """
        
        return_files=[("xtb.trj", out_file)]
        md = xtb_run(self.xtb_bin, geom_file,
                     "--metadyn",
                     *self.extra_args,
                     xcontrol=xcontrol,
                     prefix="MTD",
                     delete=self.delete,
                     scratch=self.scratchdir,
                     failout=failout,
                     return_files=return_files)
        return md




