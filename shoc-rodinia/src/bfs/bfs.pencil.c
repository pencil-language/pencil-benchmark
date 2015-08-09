#include <pencil.h>

void parallel_version( int no_of_nodes, int edge_start_no[const restrict static no_of_nodes], int edge_count[const restrict static no_of_nodes],
                       int no_of_edges, int dst_node_index[ const restrict static no_of_edges], int cost[const restrict static no_of_edges],
                       char front[const restrict static no_of_nodes], char updating_front[const restrict static no_of_nodes],
                       char visited[const restrict static no_of_nodes], int src_index)
{
  unsigned int carry_on = 1;

  for(int i = 0; i < no_of_nodes; i++ ) {
    front[i]          = 0 ;
    updating_front[i] = 0;
    visited[i]        = 0;
    cost[i]           = -1;
  }

  front[   src_index ] = 1;
  visited[ src_index ] = 1;
  cost[    src_index ] = 0;


  while( carry_on == 1 ) {
    carry_on = 0;
#pragma scop
    {
      __pencil_assume (no_of_nodes >= 2);
      __pencil_assume (no_of_edges >= 1);
#pragma pencil independent
      for( int i = 0; i < no_of_nodes; i++){
        if( front[i] == 1) {
          front[i] = 0;
          for( int j = edge_start_no[i]; j < (edge_start_no[i] + edge_count[i]) ; j++ ){
            int dst_node = dst_node_index[j];
            if( visited[ dst_node ] == 0 )
            {
              cost[dst_node] = cost[ i ] + 1;

              updating_front[ dst_node ] = 1;
            }
          }
        }
      }

#pragma pencil independent
      for( int i= 0; i < no_of_nodes; i++ ) {
        if( updating_front[i] == 1 ){
          front[i]     = 1;
          visited[i]  = 1;
          carry_on  = 1;
          updating_front[i] = 0 ;
        }
      }
    }
#pragma endscop
  }
}
