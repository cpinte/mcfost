#include "voro++.hh"
//#include <iostream>
using namespace voro;

extern "C" {
  void voro_C(int n_points, double limits[6], double x[n_points], double y[n_points], double z[n_points],
            int &n_in, double volume[], int first_neighbours[], int last_neighbours[], int &n_neighbours_tot, int neighbours_list[], int &ierr) {

    //std::cout << "Hello World!" << n_points << std::endl;

    ierr = 0 ;

    // intent(in)
    double ax,bx, ay,by, az, bz ;
    ax = limits[0] ; bx = limits[1] ;
    ay = limits[2] ; by = limits[3] ;
    az = limits[4] ; bz = limits[5] ;

    int i, nx,ny,nz,init_mem(8);
    //wall_list wl; // I am adding extra walls yet
    bool xperiodic,yperiodic,zperiodic ;
    pre_container *pcon=NULL; // pre-container pointer

    xperiodic = false ;
    yperiodic = false ;
    zperiodic = false ;

    // Define a pre-container to determine the optimal block size
    pcon=new pre_container(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
    for(i=0;i<n_points;i++) {
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
    double xx,yy,zz,r;
    std::vector<int> vi;   // vi.size() pour avoir le nombre d'elements

    int n_neighbours, first_neighbour, last_neighbour ;
    int max_size_list = 20 * n_points ; // todo : manage to get the size of the fortran array

    n_neighbours_tot = 0 ;
    last_neighbour = -1 ;
    n_in = 0 ;
    if(vlo.start()) do if(con.compute_cell(c,vlo)) { // return fals if the cell was removed
          //std::cout << "test1" << std::endl;
          n_in++ ;
          //std::cout << n_in << std::endl;

          // id of the current cell in the c_loop
          pid = vlo.pid() ;

          // Volume
          //std::cout << pid << std::endl;

          //std::cout << "test2" << std::endl;
          volume[pid] = c.volume() ;

          // Store the neighbours list
          n_neighbours = c.number_of_faces() ;
          n_neighbours_tot = n_neighbours_tot + n_neighbours ;

          first_neighbour = last_neighbour+1 ; first_neighbours[pid] = first_neighbour ;
          last_neighbour  = last_neighbour + n_neighbours ; last_neighbours[pid]  = last_neighbour ;

          //std::cout << "test3" << std::endl;

          //std::cout << n_neighbours_tot << " AA "<< max_size_list << std::endl;
          if (n_neighbours_tot > max_size_list) {
            ierr = 1 ;
            exit(1) ;
          }

          //std::cout << "test4" << std::endl;

          c.neighbors(vi) ;
          for (i=0 ; i<n_neighbours ; i++) {
           // std::cout << "**************************************" << std::endl;
           // std::cout << "A " << i << " " << vi[i] << std::endl;
           // std::cout << "fn = " << first_neighbour << std::endl;
           // std::cout << "B " << i << " " << neighbours_list[first_neighbour+1+i] << std::endl;
            neighbours_list[first_neighbour+i] = vi[i] ;
          }
          //std::cout << "test5" << std::endl;


        } while(vlo.inc());
    //std::cout << "END C! " << ierr << " " << n_in << std::endl;
  }
}
