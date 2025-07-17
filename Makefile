parallel: main_parallel.f90
	caf -cpp -DPARALLEL main_parallel.f90 -o shallow_parallel.exe

serial: main_parallel.f90
	gfortran -cpp -DSERIAL main_parallel.f90 -o shallow_serial.exe

run_serial:
	./shallow_serial.exe

run_parallel:
	cafrun -n $(CORES) ./shallow_parallel.exe
