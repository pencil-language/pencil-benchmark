#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "prl.h"
#include "measure-time.h"

//in pencil.c
void parallel_version( int no_of_nodes, int edge_start_no[no_of_nodes], int edge_count[no_of_nodes], int no_of_edges,
    int dst_node_index[ no_of_edges], int cost[no_of_edges], char front[no_of_nodes], char updating_front[no_of_nodes],
    char visited[no_of_nodes], int src_index);

void serial_version( int no_of_nodes, int edge_start_no[no_of_nodes], int edge_count[no_of_nodes], int no_of_edges,
    int dst_node_index[ no_of_edges], int cost[no_of_edges], char front[no_of_nodes], char updating_front[no_of_nodes],
    char visited[no_of_nodes], int src_index)
{

  for(int i = 0; i < no_of_nodes; i++ ) {
    front[i]          = 0;
    updating_front[i] = 0;
    visited[i]        = 0;
    cost[i]           = -1;
  }

  front[   src_index ] = 1;
  visited[ src_index ] = 1;
  cost[    src_index ] = 0;

  unsigned int carry_on = 1;
  while(carry_on) {
    carry_on = 0;

    for( int i = 0; i < no_of_nodes; i++){
      if(front[i] == 1){
        front[i] = 0;
        for( int j = edge_start_no[i]; j < (edge_start_no[i] + edge_count[i]) ; j++ ){
          int dst_node = dst_node_index[j];
          if( !visited[ dst_node ] )
          {
            cost[dst_node] = cost[ i ] + 1;

            updating_front[ dst_node ] = 1;
          }
        }
      }
    }

    for( int i= 0; i < no_of_nodes; i++ ){
      if( updating_front[i] == 1 ){
        front[i]     = 1;
        visited[i]  = 1;
        carry_on  = 1;
        updating_front[i] = 0;
      }
    }
  }
}

/*read graph file and store in malloced arrays*/
int main( int argc, char *argv[] )
{
  FILE *fp;
  int no_of_nodes;
  int no_of_edges;
  int src_index;
  int dummy_cost;

  int check = 0;

  int ret = 0;

  if(argc < 2 ){
    fprintf(stderr, "Usage: %s graph_file_name [--check]\n", argv[0]); exit(0);
  }

  fp = fopen(argv[1],"r");

  if(!fp){
    printf("unable to open graph file - %s \n",argv[1]); exit(0);
  }

  if(argc > 2 && !strcmp(argv[2], "--check"))
  {
    check = 1;
  }

  fscanf(fp,"%d",&no_of_nodes);

  int *edge_start_no = (int *)prl_alloc( sizeof(int) * no_of_nodes); //index to edge-array. starting point.
  int *edge_count = (int *)prl_alloc(sizeof(int)*no_of_nodes);        //no of edges, vertex i has.

  for(int i = 0; i < no_of_nodes; i++ ){
    fscanf(fp, "%d %d", &edge_start_no[i], &edge_count[i]); //read vertex info
  }

  fscanf(fp,"%d",&src_index);
  src_index = 0;

  fscanf(fp,"%d", &no_of_edges);

  int *dst_node_index = (int *)prl_alloc(sizeof(int) * no_of_edges);
  int *cost           = (int *)prl_alloc(sizeof(int) * no_of_edges);
  int *cost_ref       = (int *)malloc(sizeof(int) * no_of_edges);
  char *front          = (char *)prl_alloc(sizeof(char) * no_of_nodes);
  char *updating_front = (char *)prl_alloc(sizeof(char) * no_of_nodes);
  char *visited        = (char *)prl_alloc(sizeof(char) * no_of_nodes);

  for( int i = 0; i < no_of_edges; i++ ){
    fscanf(fp,"%d", &dst_node_index[i]);                         //index of destionation node
    fscanf(fp,"%d", &dummy_cost);
  }
  if(fp) fclose(fp);

  serial_version( no_of_nodes, edge_start_no, edge_count, no_of_edges, dst_node_index, cost_ref, front, updating_front, visited, src_index);

  START_MEASURE_TIME
  parallel_version( no_of_nodes, edge_start_no, edge_count, no_of_edges, dst_node_index, cost, front, updating_front, visited, src_index);
  STOP_MEASURE_TIME

  if (check)
  {
    serial_version( no_of_nodes, edge_start_no, edge_count, no_of_edges, dst_node_index, cost_ref, front, updating_front, visited, src_index);
    int mis_match = 0;
    for( int i = 0; i < no_of_nodes; i++ ) if( cost[i] != cost_ref[i]) mis_match++;
    if( mis_match == 0 )
      printf("verified: all correct :-) \n");
    else
    {
      printf("verified: INCORRECT %d mis-matches :-( \n", mis_match);
      ret = 1;
    }
  }

  PRINT_TIME

  prl_free(edge_start_no);
  prl_free(edge_count);
  prl_free(dst_node_index);
  prl_free(cost);
  free(cost_ref);
  prl_free(front);
  prl_free(updating_front);
  prl_free(visited);
  
  return ret;
}
