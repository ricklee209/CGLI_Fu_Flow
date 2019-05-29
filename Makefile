OBJS = main.o Jacobian.o Jacobian_u.o Jacobian_v.o Jacobian_w.o Initial_condition.o Driven_force.o X_boundary_condition.o Y_boundary_condition.o Z_boundary_condition.o Flux_X.o Flux_Y.o Flux_Z.o Viscous_terms.o Dual_time_stepping.o DPLUSGS.o Output.o Output_plot3d.o Statistic.o Hyperplane.o Grid.o Filter.o

HEAD = Resolution.h ijk.h prm.h MPI_prm.h Mpi_division.h


.SUFFIXES: .o .cpp

.cpp.o:
	mpicxx -O3 -c $< 



a.out : $(OBJS)

	mpicxx -O3 -o a.out $(OBJS) 





main.o: main.h Array.h Delete_Array.h $(HEAD)

Jacobian.o: $(HEAD) 

Jacobian_u.o: $(HEAD) 
Jacobian_v.o: $(HEAD) 
Jacobian_w.o: $(HEAD) 

Initial_condition.o: $(HEAD) 
Driven_force.o: $(HEAD) 

X_boundary_condition.o: $(HEAD) 
Y_boundary_condition.o: $(HEAD) 
Z_boundary_condition.o: $(HEAD) 

Flux_X.o: $(HEAD) 
Flux_Y.o: $(HEAD) 
Flux_Z.o: $(HEAD) 

Viscous_terms.o: $(HEAD) Viscous_terms.h

Dual_time_stepping.o: $(HEAD)

DPLUSGS.o: $(HEAD)

Output.o: $(HEAD)

Output_plot3d.o: $(HEAD)

Statistic.o: $(HEAD) Viscous_terms.h

Hyperplane.o: $(HEAD) Viscous_terms.h

Grid.o: $(HEAD)

Filter.o: $(HEAD)

.PHONY: clean
clean:
	rm $(OBJS)


