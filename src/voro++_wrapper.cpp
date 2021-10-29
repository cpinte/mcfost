#include <voro++/voro++.hh>
#include <iostream>
#include <iomanip>
using namespace voro;

void progress_bar(float progress) {

  int barWidth = 50;

  std::cout  << " " << std::setfill(' ') << std::setw(3) << int(progress * 100.0) << "%";
  std::cout << " |";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << " "; //">";
    else std::cout << " ";
  }
  std::cout << "| \r";
  std::cout.flush();

  if (progress >= 1.0) std::cout << std::endl;
}

#define n_elements(x)  (sizeof(x) / sizeof((x)[0]))

int index_star(int icell, int n_stars, int *stars_id) {
	int k;	
	for (k=0; k<n_stars; k++) {
		if ( (icell==stars_id[k]) ? true : false ) return k;
	}
	return -1;
}


extern "C" {
  void voro_C(int n, int max_neighbours, double limits[6], double x[], double y[], double z[], double h[],  double threshold, int n_vectors, double cutting_vectors[][3], double cutting_distance_o_h, int icell_start, int icell_end, int cpu_id, int n_cpu, int n_points_per_cpu,
              int &n_in, double volume[], int first_neighbours[], int last_neighbours[], int n_neighbours[], int neighbours_list[], bool was_cell_cut[],int n_stars, int stellar_id[], double stellar_radius[], bool stellar_neighb[], int &ierr) {

    ierr = 0 ;

    double ax,bx, ay,by, az, bz;
    ax = limits[0]; bx = limits[1];
    ay = limits[2]; by = limits[3];
    az = limits[4]; bz = limits[5];

    int i, k, nx,ny,nz,init_mem(16);

    int row = n_points_per_cpu * max_neighbours;

    //wall_list wl; // I am not adding extra walls yet
    bool xperiodic,yperiodic,zperiodic;

    xperiodic = false;
    yperiodic = false;
    zperiodic = false;

    // Makes a guess at the optimal grid of blocks to use, as in voro++ pre_container
    double dx=bx-ax, dy=by-ay, dz=bz-az;
    double optimal_particles = 5.6;
    double ilscale=pow(n/(optimal_particles*dx*dy*dz),1/3.0);
    nx=int(dx*ilscale+1);
    ny=int(dy*ilscale+1);
    nz=int(dz*ilscale+1);
    
    //Local bool array to know if a cell is a star
    //an element per cell/per_cpu
    bool *is_a_star, find_at_least_one_star=false;
    int i_star = -1; //< 0 no stars!
    is_a_star = (bool *) malloc(n_points_per_cpu * sizeof(bool));

    // Define the container
    particle_order po;
    container con(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);
    for(i=0;i<n;i++) {
      con.put(po,i,x[i],y[i],z[i]);
    }

    // Perform the Voronoi tesselation
    c_loop_order vlo(con,po);

    voronoicell_neighbor c(con);
    int pid, pid_loc;
    std::vector<int> vi;

    int n_neighbours_cell, first_neighbour, last_neighbour;
    int max_size_list = max_neighbours * n;
    double cx, cy, cz, cutting_distance, delta_edge;

    n_neighbours[cpu_id] = 0;
    last_neighbour = -1;

    float progress = 0.0;
    float progress_bar_step = 0.01;
    float pb_threshold = progress_bar_step*(1.*n)/n_cpu;

    if (!vlo.start()) {
      std::cout << "Error : voro++ did not manage to initialize" << std::endl;
      std::cout << "Exiting" << std::endl;
      ierr = 1;
      exit(1);
    }

    pid_loc = -1; // pid_loc is the cell index for this core, ie starting at 0
    do {
      pid = vlo.pid(); // id of the current cell in the c_loop

      if ((pid >= icell_start) && (pid <= icell_end)) { // We only compute a fraction of the cells in a given thread
        if (con.compute_cell(c,vlo)) { // return false if the cell was removed
          pid_loc++;

          if (cpu_id == n_cpu-1) {
            if (pid_loc > pb_threshold) {
              progress  += progress_bar_step;
              pb_threshold += progress_bar_step*(1.*n)/n_cpu;
              progress_bar(progress);
            }
          }

          // Store the neighbours list
          n_neighbours_cell = c.number_of_faces();
          n_neighbours[cpu_id] = n_neighbours[cpu_id] + n_neighbours_cell;

          first_neighbour = last_neighbour+1; first_neighbours[pid] = first_neighbour;
          last_neighbour  = last_neighbour + n_neighbours_cell; last_neighbours[pid]  = last_neighbour;

          if (n_neighbours[cpu_id] > max_size_list) {ierr = 1; exit(1);}
          
          //Is this particle a star ?
          i_star = index_star(pid, n_stars, stellar_id);
   		  is_a_star[pid_loc] = (i_star > -1) ? true : false;
   		  if (is_a_star[pid_loc]) {
   		 	puts("find a star!");
   		 	printf("pid = %d, pid_loc = %d, i_star = %d\n",pid,pid_loc,i_star);
   		 	puts("");
   		 	if (!find_at_least_one_star) find_at_least_one_star=true;
   		 	stellar_neighb[pid_loc] = false;
   		  }

          c.neighbors(vi);
          for (i=0; i<n_neighbours_cell; i++) {
            if (vi[i] >=0) {
              neighbours_list[row * cpu_id + first_neighbour+i] = vi[i] + 1;
            } else { // Wall
              neighbours_list[row * cpu_id + first_neighbour+i] = vi[i];
            }
          }

          // Compute the maximum distance to a vertex
          delta_edge = sqrt(c.max_radius_squared()); // does not cost anything

          // If the Voronoi cell is elongated, we intersect it with a dodecahedron
          was_cell_cut[pid_loc] = false;
          if (delta_edge > threshold * h[pid]) {
            cutting_distance = cutting_distance_o_h * h[pid];

            // Adding the n-planes to cut the cell
            for (k=0; k<n_vectors; k++) {
              // All the planes are a distance of 0.5 x vector length away from the center of the cell --> extra factor 2
              c.plane(cutting_vectors[k][0], cutting_vectors[k][1], cutting_vectors[k][2], 2*cutting_distance);
            }
            was_cell_cut[pid_loc] = true;
          }
          // Volume of the cell (computed after eventual cut)
          volume[pid_loc] = c.volume();
        }  // con.compute_cell
      } // pid test
    } while(vlo.inc()); //Finds the next particle to test
    n_in = pid_loc+1;

    if (cpu_id == n_cpu-1) progress_bar(1.0);
    
    /* Loop again over all cells but only do operations on cells that are stars */
    if (find_at_least_one_star) {
    	puts("Find a star, cutting neigbours!");
    } else {
    	free(is_a_star);
    	return;
    }

    pid_loc = 0;
    do {
      pid = vlo.pid(); // id of the current cell in the c_loop

      if ((pid >= icell_start) && (pid <= icell_end)) { // We only compute a fraction of the cells in a given thread
        if (is_a_star[pid_loc]) { // return false if the cell is not a star
		  con.compute_cell(c,vlo);
		  i_star = index_star(pid, n_stars, stellar_id);

          c.neighbors(vi);
          for (i=0; i<c.number_of_faces(); i++) {
            if (vi[i] >=0) {
            //Neighbour is not a wall, cut it at the surface of the star and recompute volume

          	// Compute the maximum distance to a vertex
          	delta_edge = sqrt(c.max_radius_squared()); // does not cost anything

          	//if (delta_edge > threshold * h[pid]) {
            	cutting_distance = delta_edge - stellar_radius[i_star];

            	// Adding the n-planes to cut the cell
            	for (k=0; k<1; k++) {
             	 // All the planes are a distance of 0.5 x vector length away from the center of the cell --> extra factor 2
              		c.plane(x[vi[i]]-x[pid], y[vi[i]]-x[pid], z[vi[i]]-x[pid], 2*cutting_distance);
            	}
            	//Set to True because it was cut, but make sure later that it is not treated as a wall.
            	was_cell_cut[pid_loc] = true;
          	//}
          	// Volume of the cell (computed after eventual cut)
         	 volume[pid_loc] = c.volume();

            }
          }//loop over neighbours
          
          pid_loc++;
        }  // is a star
      } // pid test
    } while(vlo.inc()); //Finds the next particle to test
    // free temporary storage before leaving
    free(is_a_star);
  }
}
