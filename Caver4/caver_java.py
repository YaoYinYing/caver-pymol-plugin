# CAVER Copyright Notice
# ============================
#

"""
"THE BEERWARE LICENSE" (Revision 42):

Yinying wrote the refactored code.
As long as you retain this notice, you can do whatever you want with this stuff.
If we meet someday, and you think this stuff is worth it, you can buy me a beer
in return.
-- Yinying Yao

"""


import logging
import os
import shutil

from .utils.live_run import run_command


class PyJava:
    def __init__(self, customized_memory_heap, caverfolder, caverjar, outdirInputs, cfgnew, out_dir):
        self.java_bin = shutil.which("java")

        logging.info("\n*** Testing if Java is installed ***")
        if not self.java_bin:
            raise RuntimeError("Java is not installed. Please install Java and try again.")

        run_command([self.java_bin, "-version"], verbose=True)

        self.jar = caverjar

        logging.info("\n*** Optimizing memory allocation for Java ***")
        self.optimize_memory(customized_memory_heap)
        self.cmd = [
            self.java_bin,
            f"-Xmx{self.memory_heap_level}m",
            "-cp",
            os.path.join(caverfolder, "lib"),
            "-jar",
            caverjar,
            "-home",
            caverfolder,
            "-pdb",
            outdirInputs,
            "-conf",
            cfgnew,
            "-out",
            out_dir,
        ]
        logging.info("*** Caver will be called using command ***")
        logging.info(" ".join(['"%s"' % t if t != "java" and t[0] != "-" else t for t in self.cmd]))
        logging.info("******************************************")

    def run_caver(self):
        return run_command(self.cmd, verbose=True)

    def optimize_memory(self, customized_memory_heap):
        customized_memory_heap = int(customized_memory_heap)
        memory_allocate_options = [
            500,
            800,
            900,
            950,
            1000,
            1050,
            1100,
            1150,
            1200,
            1250,
            1300,
            1400,
            1500,
            2000,
            3000,
            4000,
            5000,
            6000,
            8000,
            10000,
            14000,
            16000,
            20000,
            32000,
            48000,
            64000,
        ]
        memory_allocate_options.append(customized_memory_heap)
        memory_allocate_options.sort()
        # sorted(values)
        self.memory_heap_level = memory_allocate_options[0]
        for heap_level in memory_allocate_options:
            if int(heap_level) <= customized_memory_heap:
                cmd = [self.java_bin, f"-Xmx{heap_level}m", "-jar", self.jar, "do_nothing"]
                heap_level_test = run_command(cmd)

                if heap_level_test.returncode == 0:
                    self.memory_heap_level = heap_level
                    logging.debug("Memory heap level: " + str(self.memory_heap_level))
        logging.info("*** Memory for Java: " + str(self.memory_heap_level) + " MB ***")
