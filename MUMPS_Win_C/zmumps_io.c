/*

   THIS FILE IS PART OF MUMPS VERSION 4.6.2
   This Version was built on Fri Apr 14 14:59:20 2006


  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
  Since this first public domain version in 1999, the developments are
  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
  INRIA.

  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

  Up-to-date copies of the MUMPS package can be obtained
  from the Web pages http://www.enseeiht.fr/apo/MUMPS/
  or http://graal.ens-lyon.fr/MUMPS


   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.


  User documentation of any code that uses this software can
  include this complete notice. You can acknowledge (using
  references [1], [2], and [3] the contribution of this package
  in any scientific publication dependent upon the use of the
  package. You shall use reasonable endeavours to notify
  the authors of the package of this publication.

   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
   Multifrontal parallel distributed symmetric and unsymmetric solvers,
   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).

   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
   A fully asynchronous multifrontal solver using distributed dynamic
   scheduling, SIAM Journal of Matrix Analysis and Applications,
   Vol 23, No 1, pp 15-41 (2001).

   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
   S. Pralet, Hybrid scheduling for the parallel solution of linear
   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).

*/
/*    $Id: zmumps_io.c,v 1.47 2006/04/06 09:40:29 aguermou Exp $    */

#include "zmumps_io_basic_extern.h"
#include "zmumps_io_err_extern.h"
#include "zmumps_io.h"

#ifndef _WIN32
#include "zmumps_io_thread_extern.h"

#endif

#ifndef _WIN32
double zmumps_time_spent_in_sync;
#endif

double read_op_vol,write_op_vol,total_vol;

int MUMPS_CALL zmumps_is_there_finished_request(int* flag, int* ierr){
  /* Checks if there is a finished request in the queue of finished requests */
  /* On return flag=1 if there a finished request (flag=0 otherwise)         */
#ifndef _WIN32
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  switch(zmumps_io_flag_async){
  case IO_SYNC: 
    printf("zmumps_is_there_finished_request should not be called with strategy %d\n",zmumps_io_flag_async);
    break;
#ifndef _WIN32
  case IO_ASYNC_TH:
    *ierr=zmumps_is_there_finished_request_th(flag);
    if(*ierr<0){
      zmumps_io_prop_err_info(*ierr);
    }
    break;
#endif
  default: 
    sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
    *ierr=-92;
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
#ifndef _WIN32
  gettimeofday(&end_time,NULL);
  zmumps_time_spent_in_sync=zmumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return 0;
}

int MUMPS_CALL zmumps_clean_request(int* request_id,int *ierr){
#ifndef _WIN32
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  switch(zmumps_io_flag_async){
  case IO_SYNC: 
    printf("zmumps_clean_request should not be called with strategy %d\n",zmumps_io_flag_async);
    break;
#ifndef _WIN32
  case IO_ASYNC_TH:
    *ierr=zmumps_clean_request_th(request_id);
    if(*ierr<0){
      zmumps_io_prop_err_info(*ierr);
    }
    break;
#endif
  default: 
    sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
    *ierr=-92;
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
#ifndef _WIN32
  gettimeofday(&end_time,NULL);
  zmumps_time_spent_in_sync=zmumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return 0;
}

int MUMPS_CALL zmumps_test_request(int* request_id,int *flag,int* ierr){
#ifndef _WIN32
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  switch(zmumps_io_flag_async){
  case IO_SYNC: 
    printf("zmumps_test_request should not be called with strategy %d\n",zmumps_io_flag_async);
    break;
#ifndef _WIN32
  case IO_ASYNC_TH:
    *ierr=zmumps_test_request_th(request_id,flag);
    if(*ierr<0){
      zmumps_io_prop_err_info(*ierr);
    }
    break;
#endif
  default: 
    sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
    *ierr=-92;
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }  
#ifndef _WIN32
  gettimeofday(&end_time,NULL);
  zmumps_time_spent_in_sync=zmumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return 0;
}

int MUMPS_CALL zmumps_wait_request(int *request_id,int* ierr){
#ifndef _WIN32
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  if(*request_id==-1)
    return 0;
  switch(zmumps_io_flag_async){
  case IO_SYNC: 
    printf("zmumps_wait_request should not be called with strategy %d\n",zmumps_io_flag_async);
    break;
#ifndef _WIN32
  case IO_ASYNC_TH:
    *ierr=zmumps_wait_request_th(request_id);
    if(*ierr<0){
      zmumps_io_prop_err_info(*ierr);
    }
    break;
#endif
  default: 
    sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
    *ierr=-92;
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
    /*    printf("Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
          exit (-3);*/
  }
#ifndef _WIN32
  gettimeofday(&end_time,NULL);
  zmumps_time_spent_in_sync=zmumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return 0;
}

int MUMPS_CALL zmumps_wait_all_requests(int *ierr){
#ifndef _WIN32
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  switch(zmumps_io_flag_async){
  case IO_SYNC: 
    printf("zmumps_wait_all_requests should not be called with strategy %d\n",zmumps_io_flag_async);
    break;
#ifndef _WIN32
  case IO_ASYNC_TH:
    *ierr=zmumps_wait_all_requests_th();
    if(*ierr<0){
      zmumps_io_prop_err_info(*ierr);
    }
    break;
#endif
  default: 
    sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
    *ierr=-92;
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
#ifndef _WIN32
  gettimeofday(&end_time,NULL);
  zmumps_time_spent_in_sync=zmumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return 0;
}

int MUMPS_CALL zmumps_low_level_init_ooc_c(int* _myid, int* total_size_io,int* size_element,
                                          int* async, int* k211, char* zmumps_dir, char* zmumps_file,
                                          int* zmumps_dim_dir, int* zmumps_dim_file,
                                          int* ierr, zmumps_ftnlen l1, zmumps_ftnlen l2){
  /* Computes the number of files needed. Uses ceil value. */
  zmumps_io_nb_file=0;
  zmumps_io_last_file_opened=-1;
#ifdef _WIN32
  printf("inside %d\n",*async);
  if(*async==IO_ASYNC_AIO||*async==IO_ASYNC_TH){
    zmumps_io_is_init_called=0;
    sprintf(error_str,"Error: Forbidden value of Async flag with _WIN32\n");
    *ierr=-92;
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
#endif
  total_vol=0;
  zmumps_io_flag_async=*async;
  zmumps_io_k211=*k211;
  *ierr=zmumps_init_file_name(zmumps_dir,zmumps_file,zmumps_dim_dir,zmumps_dim_file,_myid);
  if(*ierr<0){
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
  *ierr=zmumps_init_file_structure(_myid,total_size_io,size_element);
  if(*ierr<0){
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
#ifndef _WIN32
  zmumps_time_spent_in_sync=0;
#endif


  if(*async){
    switch(*async){
    case IO_SYNC: 
      printf("zmumps_low_level_init_ooc_c should not be called with strategy %d\n",zmumps_io_flag_async);
      break;
#ifndef _WIN32
    case IO_ASYNC_TH:
      zmumps_low_level_init_ooc_c_th(async,ierr);
      if(*ierr<0){
        zmumps_io_prop_err_info(*ierr);
        return *ierr;
      }
      break;
#endif
    default: 
      /*      printf("Error: unknown I/O strategy : %d\n",*async);*/
      *ierr=-92;
      zmumps_io_prop_err_info(*ierr);
      return *ierr;
      /*      exit (-3);*/
    }
  }
  zmumps_io_is_init_called=1;
  return 0;
}


/**
 * Writes a contigous block of central memory to the disk.
 */
int MUMPS_CALL zmumps_low_level_write_ooc_c(const int * strat_IO, 
                                void * address_block,
                                int * block_size,
                                int * pos_in_file,
                                int * file_number,
                                int * inode,
                                int * request_arg,
                                int * ierr){
   int ret_code=0;
#ifndef _WIN32
   struct timeval start_time,end_time;
   gettimeofday(&start_time,NULL);
#endif
   if(zmumps_io_flag_async){
     switch(*strat_IO){
#ifndef _WIN32
     case IO_ASYNC_TH:
       ret_code=zmumps_async_write_th(strat_IO, address_block, block_size,pos_in_file,file_number,inode,request_arg,ierr);       
       if(ret_code<0){
         *ierr=ret_code;
         zmumps_io_prop_err_info(ret_code);
       }
       break;
#endif
     default: 
       ret_code=-91;
       *ierr=ret_code;
       sprintf(error_str,"Error: unknown I/O strategy : %d\n",*strat_IO);
       zmumps_io_prop_err_info(ret_code);
       return ret_code;
     }
   }else{
     ret_code=zmumps_io_do_write_block(address_block,block_size,pos_in_file,file_number,ierr);   
     if(ret_code<0){
       *ierr=ret_code;
       zmumps_io_prop_err_info(ret_code);
     }
   }
#ifndef _WIN32
   gettimeofday(&end_time,NULL);
   zmumps_time_spent_in_sync=zmumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
   write_op_vol=write_op_vol+((*block_size)*zmumps_elementary_data_size);
   return ret_code;
}

/**
 * Writes a contigous block of central memory to the disk.
 **/
int MUMPS_CALL zmumps_low_level_read_ooc_c(const int * strat_IO, 
                                          void * address_block,
                                          int * block_size,
                                          int * from_where,
                                          int * file_number,
                                          int * inode,
                                          int * request_arg,
                                          int * ierr){
  int ret_code=0;
#ifndef _WIN32
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  if(zmumps_io_flag_async){
      switch(*strat_IO){
#ifndef _WIN32
      case IO_ASYNC_TH:
        ret_code=zmumps_async_read_th(strat_IO,address_block,block_size,from_where,file_number,inode,request_arg,ierr);
        if(ret_code<0){
          *ierr=ret_code;
          zmumps_io_prop_err_info(ret_code);
        }
        break;
#endif
      default: 
        ret_code=-91;
        *ierr=ret_code;
        sprintf(error_str,"Error: unknown I/O strategy : %d\n",*strat_IO);
        zmumps_io_prop_err_info(ret_code);
        return ret_code;
      }
  }else{
    ret_code=zmumps_io_do_read_block(address_block,block_size,from_where,file_number,ierr);
    if(ret_code<0){
      *ierr=ret_code;
      zmumps_io_prop_err_info(ret_code);
    }
    *request_arg=1;
  }
#ifndef _WIN32
  gettimeofday(&end_time,NULL);
  zmumps_time_spent_in_sync=zmumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  read_op_vol=read_op_vol+((*block_size)*zmumps_elementary_data_size);
  return ret_code;
}

int MUMPS_CALL zmumps_low_level_direct_read(void * address_block,
                                int * block_size,
                                int * from_where,
                                int * file_number,
                                int * ierr){
  int ret_code=0;
  
#ifndef _WIN32
    if(zmumps_io_flag_async==IO_ASYNC_TH||zmumps_io_flag_async==0){
#else
    if(zmumps_io_flag_async==0){
#endif
      ret_code=zmumps_io_do_read_block(address_block,block_size,from_where,file_number,ierr);
      if(ret_code<0){
         *ierr=ret_code;
	 zmumps_io_prop_err_info(ret_code);
	 return ret_code;
      }
    }else{
    }
    read_op_vol=read_op_vol+((*block_size)*zmumps_elementary_data_size);
    return ret_code;
}

  
int MUMPS_CALL zmumps_clean_io_data_c(int * myid,int *ierr){
  /* cleans the thread/io management data*/
  if(!zmumps_io_is_init_called) return -1;
  switch(zmumps_io_flag_async){
  case IO_SYNC: 
    break;
#ifndef _WIN32
  case IO_ASYNC_TH:
    *ierr=zmumps_clean_io_data_c_th(myid);
    if(*ierr<0){
      zmumps_io_prop_err_info(*ierr);
    }
    break;
#endif
  default: 
    *ierr=-91;
    sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }

  zmumps_free_file_pointers();
  zmumps_io_is_init_called=0;

  return 0;
}

int MUMPS_CALL zmumps_ooc_print_stats(){
#ifndef _WIN32
  printf("%d: total time spent in i/o mode = %lf\n",zmumps_io_myid,zmumps_time_spent_in_sync);
#endif
  printf("%d: Volume of read i/o = %lf\n",zmumps_io_myid,read_op_vol);
  printf("%d: Volume of write i/o = %lf\n",zmumps_io_myid,write_op_vol);
  total_vol=total_vol+read_op_vol+write_op_vol;
  printf("%d: Total i/o volume = %lf\n",zmumps_io_myid,total_vol);
  return 0; 
}
int MUMPS_CALL zmumps_get_max_nb_req(int *max,int *ierr){
  switch(zmumps_io_flag_async){
  case IO_SYNC: 
    *max=1;
    break;
#ifndef _WIN32
  case IO_ASYNC_TH:
    *max=MAX_FINISH_REQ+MAX_IO;
    break;
#endif
  default: 
    *ierr=-91;
    sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
    zmumps_io_prop_err_info(*ierr);    
    return *ierr;
  }
  return 0;
}

int MUMPS_CALL zmumps_get_max_file_size(double * max_ooc_file_size){
  *max_ooc_file_size=(double)(MAX_FILE_SIZE);
  return 0;
}

int MUMPS_CALL zmumps_ooc_get_nb_files(int* nb_files){
  return zmumps_io_get_nb_files(nb_files);
}

int MUMPS_CALL zmumps_ooc_get_file_name(int* indice,char* name,int* length, zmumps_ftnlen l1){
  return zmumps_io_get_file_name(indice,name,length);
}

int MUMPS_CALL zmumps_ooc_set_file_name(int* indice,char* name,int* length,int *ierr, zmumps_ftnlen l1){
  *ierr=zmumps_io_set_file_name(indice,name,length);
  if(*ierr<0){
    zmumps_io_prop_err_info(*ierr);
  }
  return *ierr;
}

int MUMPS_CALL zmumps_ooc_alloc_pointers(int* dim,int* ierr){
  int ret_code;
  ret_code=zmumps_io_alloc_pointers(dim);
  if(ret_code<0){
    *ierr=ret_code;
    zmumps_io_prop_err_info(ret_code);
  }
  *ierr=ret_code;
  zmumps_io_set_last_file(dim);
  return ret_code;
}


int MUMPS_CALL zmumps_ooc_init_vars(int* myid_arg, int* nb_file_arg,
                        int* size_element,int* async, int* k211,
                        char* zmumps_dir, char* zmumps_file,
                        int* zmumps_dim_dir, int* zmumps_dim_file,
                        int *ierr, zmumps_ftnlen l1, zmumps_ftnlen l2){
#ifndef _WIN32
  zmumps_time_spent_in_sync=0;
#endif
  zmumps_io_k211=*k211;
  *ierr=zmumps_init_file_name(zmumps_dir,zmumps_file,zmumps_dim_dir,zmumps_dim_file,myid_arg);
  if(*ierr<0){
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
  return zmumps_io_init_vars(myid_arg, nb_file_arg,size_element,async);
}

int MUMPS_CALL zmumps_ooc_start_low_level(int *ierr){

  read_op_vol=0;
  write_op_vol=0;
  *ierr=zmumps_io_open_files_for_read();
  if(*ierr<0){ 
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
  if(zmumps_io_flag_async){
    switch(zmumps_io_flag_async){
    case IO_SYNC: 
      break;
#ifndef _WIN32
    case IO_ASYNC_TH:
      zmumps_low_level_init_ooc_c_th(&zmumps_io_flag_async,ierr);
      if(*ierr<0){
        zmumps_io_prop_err_info(*ierr);
        return *ierr;
      }
      break;
#endif
    default: 
      *ierr=-91;
      sprintf(error_str,"Error: unknown I/O strategy : %d\n",zmumps_io_flag_async);
      zmumps_io_prop_err_info(*ierr);    
      return *ierr;
    }
  }
  zmumps_io_is_init_called=1;
  return *ierr;
}

int MUMPS_CALL zmumps_ooc_remove_file(char *name,int *ierr, zmumps_ftnlen l1){
  *ierr=remove(name);
  if(*ierr<0){
#ifndef _WIN32
    zmumps_io_build_err_str(errno,-90,"Unable to remove OOC file",error_str,200);
#else
    sprintf(error_str,"Unable to remove OOC file");
#endif
    *ierr=-90;
    zmumps_io_prop_err_info(*ierr);
    return *ierr;
  }
  return 0;
}

int MUMPS_CALL zmumps_ooc_end_write(int *ierr){
  return 0;
}
