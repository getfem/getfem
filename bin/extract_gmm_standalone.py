#!/usr/bin/env python3
# -*- python -*-
#
# Copyright (C) 2026 Yves Renard, Konstantinos Poulios.
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program. If not, see https://www.gnu.org/licenses/.
#
############################################################################

import os
import shutil
import subprocess
from pathlib import Path

def main():
  root = Path.cwd() / "../gmm_standalone_temp"
  getfem_root = Path(__file__).parent.parent

  # Create directory structure
  shutil.rmtree(root, ignore_errors=True)
  root.mkdir()
  (root / "tests").mkdir()
  (root / "include").mkdir()
  (root / "include" / "gmm").mkdir()

  # Copy gmm header files
  for h_file in (getfem_root / "src" / "gmm").glob("gmm*.h"):
    shutil.copy2(h_file, root / "include" / "gmm" / h_file.name)
  shutil.copy2(getfem_root / "src" / "gmm" / "gmm_arch_config.h.in",
               root / "include" / "gmm" / "gmm_arch_config.h.in")
  # Copy test files
  for cc_file in (getfem_root / "tests").glob("gmm_torture*.cc"):
    shutil.copy2(cc_file, root / "tests" / cc_file.name)
  shutil.copy2(getfem_root / "tests" / "make_gmm_test.pl",
               root / "tests" / "make_gmm_test.pl")

  # Copy other files
  shutil.copy2(getfem_root / "autogen.sh",       root / "autogen.sh")
  shutil.copy2(getfem_root / "gmm-config.in",    root / "gmm-config.in")
  shutil.copy2(getfem_root / "COPYING",          root / "COPYING")
  shutil.copy2(getfem_root / "README",           root / "README")
  shutil.copy2(getfem_root / "NEWS",             root / "NEWS")
  shutil.copy2(getfem_root / "ChangeLog",        root / "ChangeLog")
  shutil.copy2(getfem_root / "configure_gmm.ac", root / "configure.ac") # renaming
  shutil.copy2(getfem_root / "Makefile_gmm.am",  root / "Makefile.am")
  shutil.copy2(getfem_root / "src" / "Makefile_gmm.am",  root / "include" / "Makefile.am")
  shutil.copy2(getfem_root / "tests" / "Makefile_gmm.am",  root / "tests" / "Makefile.am")
  # Copy m4 directory
  shutil.copytree(getfem_root / "m4", root / "m4")

  # Create AUTHORS file and dummy test file
  os.chdir(getfem_root)
  with open(root / "AUTHORS", "w") as f:
    f.write("Authors of GMM\n"
            "Yves RENARD. Initial project. All the project.\n"
            "Julien POMMIER. All the project.\n"
            "Konstantinos POULIOS. All the project.\n")
  with open(root / "tests" / "dummy.cc", "w") as f:
    f.write("#include <iostream>\nint main(void) { return 0; }\n")

  # Run autogen, configure, and make dist, using subprocess
  subprocess.run(f"cd {root} && ./autogen.sh", shell=True, check=True)
  subprocess.run(f"cd {root} && ./configure", shell=True, check=True)
  subprocess.run(f"cd {root} && make dist", shell=True, check=True)
  # Move tar.gz file
  
  tar_files = list(getfem_root.glob("gmm-*.tar.gz"))
  for tar_file in tar_files:
    os.remove(str(tar_file))

  tar_files = list(root.glob("gmm-*.tar.gz"))
  for tar_file in tar_files:
    shutil.move(str(tar_file), getfem_root)

  # Clean up
  shutil.rmtree(root)

if __name__ == "__main__":
  main()
