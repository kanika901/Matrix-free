# Notes for running Bratu (ex5s), Driven cavity (ex19)
# How to run: ./notes.sh > notes.log
# Then when you want to get the runtimes you can do
# grep -E '(Main Stage:)|(ex5 with)|(ex19 with)|(Failed to converge)' notes.log
# TODO: Check for convergence! e.g. cg with ex19 does not converge!

# How to build PETSc
#cd $HOME/rnet/petsc-3.8.4
#./configure --prefix=$HOME/rnet/build --CFLAGS="-O3 -g" --download-triangle=1 --with-debugging=no
#make PETSC_DIR=$PETSC_DIR PETSC_ARCH=arch-linux2-c-opt all
#make PETSC_DIR=$PETSC_DIR PETSC_ARCH=arch-linux2-c-opt install


# Don't load petsc since we use the old version
if [ "$HOSTNAME" = 'talapas-ln1' ]; then
	module load matlab
	module load intel/17
	module load mkl
	module load openmpi/2.1
	export PETSC_DIR=$HOME/rnet/petsc-3.8.4
fi

# Parameters (explained below)
EX19GRID=100
LAMBDA="5.5"
GRIDEX5=100
LIDVELOCITY=10
GRASHOF=1000

# ex19 Parameters:
# -prandtl: Doesn't make sense to change
# -grashof: try 100, 1000, 10000, larger means harder to solve
# -lidvelocity 10
# ex5 Parameter: -par: Between 0 and 6.81 (lambda, a.k.a. nonlinearity parameter)

# ex5 and ex19 are based on
# /projects/cis607hpc/shared/packages/petsc/3.8.3/src/petsc-3.8.3/src/snes/examples/tutorials
VARNAME=drivencavity_jac_${EX19GRID}x${EX19GRID}_lv${LIDVELOCITY}_grashof$GRASHOF
make ex19
MPIEXEC='' # Make sure to add a space after the command so this works if you leave it unset
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -dm_view ascii -da_grid_x $EX19GRID -da_grid_y $EX19GRID -ksp_monitor_jacobian
# Generate jacobian
mv ex19_mat_snes0000_ksp0000.m "$VARNAME".m
matlab -nojvm -nodisplay -nosplash -r "run('$VARNAME.m'); $VARNAME = spconvert(zzz); save('$VARNAME.mat', '$VARNAME'); exit;"
sed -i 's/Vec_.* = \[/Vec_ex19_0 = \[/' ex19_x0.m

# Run the timing experiments

# Bratu
# MMS: Method of Manufacturing Solutions
# http://www.math.ttu.edu/~klong/Sundance/html/NonlinearExamples.pdf
# Don't worry too much about m_par and n_par
VARNAME=bratu_${GRIDEX5}x${GRIDEX5}_lambda${LAMBDA/\./pt}
make ex5
${MPIEXEC}./ex5 -snes_monitor -par $LAMBDA -da_grid_x $GRIDEX5 -da_grid_y $GRIDEX5 -ksp_monitor_jacobian
mv ex5_mat_snes0000_ksp0000.m "$VARNAME".m
matlab -nojvm -nodisplay -nosplash -r "run('$VARNAME.m'); $VARNAME = spconvert(zzz); save('$VARNAME.mat', '$VARNAME'); exit;"
sed -i 's/Vec_.* = \[/Vec_ex5_0 = \[/' ex5_x0.m

# Note: ASM(ovl) means set -pc_asm_overlap <ovl> (1 is the default)
# For bratu_400x400_lambda2pt5.mat with the ML model, we get 5 good solvers-PC combinations:
# 1. LSQR + Block Jacobi 
# 2. FGMRES + ASM(1) 
# 3. iBCGS + ASM(1) 
# 4. iBCGS + ASM(0) 
# 5. iBCGS + Block Jacobi

# For the bratu_100x100_lambda5pt5 below are the ML model suggestions: 
# 1. iBCGS + ASM(1)
# 2. iBCGS + ASM(0)
# 3. iBCGS + Block Jacobi

# Ex19 100x100
#lgmres,bjacobi 
#cg,asm(0)
#lsqr,bjacobi
#gmres,asm(0) 
#fgmres,asm(1) 
#ibcgs,asm(1)
#ibcgs,asm(0)
#ibcgs,bjacobi
mkdir -p logs
echo "ex19 with gmres bjacobi"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type gmres -pc_type bjacobi > logs/ex19gmresbjacobi.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with lgmres bjacobi"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type lgmres -pc_type bjacobi > logs/ex19lgmresbjacobi.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with cg asm(0). ****NOTE: This fails to converge! ****"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type cg -pc_type asm -pc_asm_overlap 0 > logs/ex19cgasm0.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with lsqr bjacobi"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type lsqr -pc_type bjacobi > logs/ex19lsqrbjacobi.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with gmres asm(0)"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type gmres -pc_type asm -pc_asm_overlap 0 > logs/ex19gmresasm0.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with fgmres asm(1)"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type fgmres -pc_type asm -pc_asm_overlap 1 > logs/ex19fgmresasm1.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with ibcgs asm(1)"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type ibcgs -pc_type asm -pc_asm_overlap 1 > logs/ex19ibcgsasm1.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with ibcgs asm(0)"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type ibcgs -pc_type asm -pc_asm_overlap 0 > logs/ex19ibcgsasm0.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

echo "ex19 with ibcgs bjacobi"
${MPIEXEC}./ex19 -snes_monitor -lidvelocity $LIDVELOCITY -prandtl 1.0 -grashof $GRASHOF -da_grid_x $EX19GRID -da_grid_y $EX19GRID -log_view  -ksp_type ibcgs -pc_type bjacobi > logs/ex19ibcgsbjacobi.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi


# Example 5
echo "ex5 with gmres bjacobi"
${MPIEXEC}./ex5 -snes_monitor -par $LAMBDA -da_grid_x $GRIDEX5 -da_grid_y $GRIDEX5 -log_view  -ksp_type gmres -pc_type bjacobi > logs/ex5gmresbjacobi.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

#echo "ex5 with lsqr bjacobi"
#${MPIEXEC}./ex5 -snes_monitor -par $LAMBDA -da_grid_x $GRIDEX5 -da_grid_y $GRIDEX5 -log_view  -ksp_type lsqr -pc_type bjacobi
#if [ "$?" -ne 0 ]; then
#	echo "Failed to converge"
#fi

#echo "ex5 with fgmres asm(1)"
#${MPIEXEC}./ex5 -snes_monitor -par $LAMBDA -da_grid_x $GRIDEX5 -da_grid_y $GRIDEX5 -log_view  -ksp_type fgmres -pc_type asm -pc_asm_overlap 1
#if [ "$?" -ne 0 ]; then
#	echo "Failed to converge"
#fi

echo "ex5 with ibcgs asm(1)"
${MPIEXEC}./ex5 -snes_monitor -par $LAMBDA -da_grid_x $GRIDEX5 -da_grid_y $GRIDEX5 -log_view  -ksp_type ibcgs -pc_type asm -pc_asm_overlap 1 > logs/ex5ibcgsasm1.log
if [ "$?" -ne 0 ]; then
	echo "Failed to converge"
fi

#echo "ex5 with ibcgs asm(0)"
#${MPIEXEC}./ex5 -snes_monitor -par $LAMBDA -da_grid_x $GRIDEX5 -da_grid_y $GRIDEX5 -log_view  -ksp_type ibcgs -pc_type asm -pc_asm_overlap 0
#if [ "$?" -ne 0 ]; then
#	echo "Failed to converge"
#fi

#echo "ex5 with ibcgs bjacobi"
#${MPIEXEC}./ex5 -snes_monitor -par $LAMBDA -da_grid_x $GRIDEX5 -da_grid_y $GRIDEX5 -log_view  -ksp_type ibcgs -pc_type bjacobi
#if [ "$?" -ne 0 ]; then
#	echo "Failed to converge"
#fi

