# 2D-Steady state diffusion equation using  Nodal Integral Method, which is a coarse mesh method and able to produce result with reasonably coarse grid size.

COPYRIGHT

This program is based on Nodal Integral Method for the 2-Dimesional Steady State Diffusion equation. 
Copyright (C) 2017 Rishabh Prakash Sharma. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
\***************************************************************\
 
HOW_TO_USE

FOR_BUILD_&_RUN_USE_FOLLOWING_COMMANDS
\*******USE_GCC_FOR_COMPILATION*******\

1. make         \**(avoid warnings)**\
2. make run     \**(run the build file)**\


FOR_PLOTTING_RESULT_USE_"MATLAB"
 
\*** S-averaged results*****\ 

\***Contour****\
1. load Xs
2. load Ys
3. load Ts
4. contourf(Xs,Ys,Ts)

\***Temperature profile***\
1. load T_P1
2. r=T_P1(:,1);
3. Ts=T_P1(:,2);
4. plot(r,Ts,'-*');


FOLLOW_THE_SAME_PROCEDURE_FOR_(T-averaged_results)

DELETE_Generated_files
1. make clean
