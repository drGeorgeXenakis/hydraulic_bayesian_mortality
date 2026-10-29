FC = gfortran
FFLAGS = -std=f95 -Wall -Wextra -O2 -fcheck=all

hydraulic_mortality: hydraulics.f90 mortality.f90
	$(FC) $(FFLAGS) hydraulics.f90 mortality.f90 -o hydraulic_mortality

run: hydraulic_mortality
	./hydraulic_mortality

check:
	python3 verify.py

clean:
	rm -f hydraulic_mortality *.o *.mod
