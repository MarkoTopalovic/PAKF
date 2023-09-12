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
/*    $Id: dmumps_io_basic.h,v 1.45 2006/04/13 11:54:15 jylexcel Exp $  */

#define MAX_FILE_SIZE 1879048192 /* (2^31)-1-(2^27) */


/*                                                      */
/* Important Note :                                     */
/* ================                                     */
/* On GNU systems, __USE_GNU must be defined to have    */
/* access to the O_DIRECT I/O flag.                     */
/*                                                      */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>

#endif

#ifndef _WIN32
#define DMUMPS_IO_FLAG_TRY_O_DIRECT O_RDWR
#endif



/* Force WITH_PFUNC on architectures where we know that it should work */
#if defined(IRIX64_) | defined(SP_) | defined(SUN_) | defined(_GNU_SOURCE)
#undef WITH_PFUNC
#define WITH_PFUNC
#endif


#define IO_SYNC      0

#define IO_ASYNC_TH  1

#define IO_ASYNC_AIO 2

#define IO_READ 1
#define IO_WRITE 0

#define UNITIALIZED "NAME_NOT_INITIALIZED"

#define DMUMPS_OOC_DEFAULT_DIR "/tmp"

#define SEPARATOR "/"

int _dmumps_next_file();
void _dmumps_update_current_file_position();
int _dmumps_compute_where_to_write(const double to_be_written);
int _dmumps_prepare_pointers_for_write(double to_be_written,int * pos_in_file, int * file_number);
int dmumps_io_do_write_block(void * address_block,int * block_size,int * pos_in_file,int * file_number,int * ierr);
int dmumps_io_do_read_block(void * address_block,int * block_size,int * from_where,int * file_number,int * ierr);
int _dmumps_compute_nb_concerned_files(int * block_size, int * nb_concerned_files);
int dmumps_free_file_pointers();
int dmumps_init_file_structure(int* _myid, int* total_size_io,int* size_element);
int dmumps_init_file_name(char* dmumps_dir,char* dmumps_file,int* dmumps_dim_dir,int* dmumps_dim_file,int* _myid);
int dmumps_io_alloc_file_struct(int* nb);
int dmumps_io_get_nb_files(int* nb_files);
int dmumps_io_get_file_name(int* indice,char* name,int* length);
int dmumps_io_alloc_pointers(int * dim);
int dmumps_io_init_vars(int* myid_arg, int* nb_file_arg,int* size_element,int* async_arg);
int dmumps_io_set_file_name(int* indice,char* name,int* length);
int dmumps_io_open_files_for_read();
int dmumps_io_set_last_file(int* dim);

int dmumps_io_write__(void *file, void *loc_add, size_t write_size, int where);

#ifndef _WIN32
int dmumps_io_write_os_buff__(void *file, void *loc_add, size_t write_size, int where);
int dmumps_io_write_direct_io__(void *file, void *loc_addr, size_t write_size, int where);
int dmumps_io_flush_write__();
#else
int dmumps_io_write_win32__(void *file, void *loc_add, size_t write_size, int where);
#endif


int dmumps_io_read__(void * file,void * loc_addr,size_t size,int local_offset);
#ifndef _WIN32
int dmumps_io_read_os_buff__(void * file,void * loc_addr,size_t size,int local_offset);
int dmumps_io_read_direct_io__(void * file,void * loc_addr,size_t size,int local_offset);
#else
int dmumps_io_read_win32__(void * file,void * loc_addr,size_t size,int local_offset);
#endif



#ifndef _WIN32  
#ifdef WITH_PFUNC

int dmumps_io_protect_pointers();
int dmumps_io_unprotect_pointers();
int dmumps_io_init_pointers_lock();
int dmumps_io_destroy_pointers_lock();

#endif /*WITH_PFUNC*/
#endif /*_WIN32*/


