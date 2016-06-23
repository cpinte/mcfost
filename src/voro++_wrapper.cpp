#include "voro++.hh"
using namespace voro;

extern "C" {
  void voro_C(int n, int max_neighbours, double limits[6], double x[], double y[], double z[],
              int &n_in, double volume[], int first_neighbours[], int last_neighbours[], int &n_neighbours_tot, int neighbours_list[], int &ierr) {

    ierr = 0 ;

    // intent(in)
    double ax,bx, ay,by, az, bz ;
    ax = limits[0] ; bx = limits[1] ;
    ay = limits[2] ; by = limits[3] ;
    az = limits[4] ; bz = limits[5] ;

    int i, nx,ny,nz,init_mem(8);

    //wall_list wl; // I am not adding extra walls yet
    bool xperiodic,yperiodic,zperiodic ;
    pre_container *pcon=NULL; // pre-container pointer

    xperiodic = false ;
    yperiodic = false ;
    zperiodic = false ;

    // Define a pre-container to determine the optimal block size
    pcon=new pre_container(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
    for(i=0;i<n;i++) {
      //std::cout << i << " " << x[i] << std::endl;
      pcon->put(i,x[i],y[i],z[i]);
    }
    pcon->guess_optimal(nx,ny,nz);

    // define the proper container and point the pre-containet toward it
    particle_order vo;
    container con(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);

    //con.add_wall(wl);
    pcon->setup(vo,con);
    delete pcon;

    c_loop_order vlo(con,vo);

    // Perform the Voronoi tesselation
    voronoicell_neighbor c;
    int pid ;
    std::vector<int> vi;

    int n_neighbours, first_neighbour, last_neighbour ;
    int max_size_list = max_neighbours * n ;

    n_neighbours_tot = 0 ;
    last_neighbour = -1 ;
    n_in = 0 ;


    if(vlo.start()) do if(con.compute_cell(c,vlo)) { // return fals if the cell was removed
          n_in++ ;

          // id of the current cell in the c_loop
          pid = vlo.pid() ;

          // Volume
          volume[pid] = c.volume() ;

          // Store the neighbours list
          n_neighbours = c.number_of_faces() ;
          n_neighbours_tot = n_neighbours_tot + n_neighbours ;

          first_neighbour = last_neighbour+1 ; first_neighbours[pid] = first_neighbour ;
          last_neighbour  = last_neighbour + n_neighbours ; last_neighbours[pid]  = last_neighbour ;

          if (n_neighbours_tot > max_size_list) {
            ierr = 1 ;
            exit(1) ;
          }

          c.neighbors(vi) ;
          for (i=0 ; i<n_neighbours ; i++) {
            if (vi[i] >=0) {
              neighbours_list[first_neighbour+i] = vi[i] + 1 ;
            } else {
              // Wall
              neighbours_list[first_neighbour+i] = vi[i] ;
            }
          }


        } while(vlo.inc()); //Finds the next particle to test
  }
}
