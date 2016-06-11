#include "voro++.hh"
using namespace voro;

//int Voronoi_tesselation() {
int main() {

  // intent(in)
  double ax,bx, ay,by, az, bz ;
  double x, y, z ; // These should be arrays

  int nx,ny,nz,init_mem(8);
  //wall_list wl; // I am adding extra walls yet
  bool xperiodic,yperiodic,zperiodic ;
  pre_container *pcon=NULL; // pre-container pointer

  xperiodic = false ;
  yperiodic = false ;
  zperiodic = false ;

  // Define a pre-container to determine the optimal block size
  pcon=new pre_container(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
  //pcon->import(argv[i+6]);  // replace by pcon.put --> liste x,y,z
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

  if(vlo.start()) do if(con.compute_cell(c,vlo)) { // return fals if the cell was removed
        //vlo.pos(pid,xz,yy,zz,r);  // renvoie pid, x, y, z
        c.number_of_faces() ;
        c.neighbors(vi) ;


  } while(vlo.inc());


}
