# Repository.  

Human 3D sinoatrial node using simple Fenton Karma cell models. 2016.

# Background.  

The sinoatrial node (SAN) is the primary natural pacemaker of the heart. This is the code repo for one of the first 3D SAN models.
The article related to this code is:  
Kharche et al. https://journals.plos.org/plosone/article/comments?id=10.1371/journal.pone.0183727

# Dependencies.  

This program solves the monodomain reaction-diffusion equations using operator splitting.
The ordinary differential equations are solved using Sundials library.
The partial differential equation (diffusion) is solved using PetSc's TS solver.
The visualization is performed using VTK/ParaView.

* PetSc can be found at: https://petsc.org/release/  
* Sundials can be found at: https://computing.llnl.gov/projects/sundials  
* ParaView can be found at: https://www.kitware.com/
* MPICH.
* GNU C.
* Python (for non-interactive viz).

# Install.  

Apart from the libraries, the C codes are compiled using provided Makefiles. Edit for locations of library sources/lib files.

# Sources and data description.  

See the following repo in my github for the geomtries used in the code:  mccsssk2 /Human-Sinoatrial-Node-Anatomical-Model  

The programs consist of a cardiomyocyte RHS and the diffusion RHS. The cardiomyoctye models are provided as individual
functions modified from the CellML website. The diffusion problem is solved using PetSc functions. In addition to an iterative solver
for the explicit finite differences, PetSc provides efficient box partitioning that allows better scaling of codes.
The output is written as binary files to help with parallel i/o.

# Use.  

The programs require a geometry, location specific cell model RHS functions. The simulation is then performed using standard
MPI commands, and PBS or slurm scripting for large computers. The outputs in binary file format is converted to VTK or information
extracted during post-processing. Python scripts that show how viz can be done as a parallel job on a cluster are provided. 
The code is extensible to any geometry and any cell model type.

# Maintainer.  

This code has now been integrated into the VCPL package.

# Acknowledements.

This project was supported by the Dobrzynski laboratory in Manchester Medical School, University of Manchester, UK. 2014-present.

# Licence.

BSD 3-Clause License

Copyright (c) 2023, Sanjay R. Kharche, Ph.D.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Dr. Sanjay R. Kharche, Ph.D. (App. Maths).  
January 23rd, 2023.  
Email: Sanjay.Kharche@lhsc.on.ca , skharche@uwo.ca , srk.london@yahoo.com .  
Phone: +15198780685 (mobile).  
Website: https://kcru.lawsonresearch.ca/research/srk/index.html  

