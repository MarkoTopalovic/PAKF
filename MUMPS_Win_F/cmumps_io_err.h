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
/*    $Id: cmumps_io_err.h,v 1.7 2006/03/14 09:49:50 jylexcel Exp $  */

#include <errno.h>

#if defined(UPPER) || defined(_WIN32)
#define cmumps_low_level_init_err_str CMUMPS_LOW_LEVEL_INIT_ERR_STR
#elif defined(Add_)
#define cmumps_low_level_init_err_str cmumps_low_level_init_err_str_
#elif defined(Add__)
#define cmumps_low_level_init_err_str cmumps_low_level_init_err_str__
#endif

#ifdef _WIN32
#define MUMPS_CALL __cdecl
#else
#define MUMPS_CALL
#endif


#ifndef _WIN32  

#include <pthread.h>

int cmumps_io_protect_err();

int cmumps_io_unprotect_err();

int cmumps_io_init_err_lock();

int cmumps_io_destroy_err_lock();

int cmumps_check_error_th();

#endif /*_WIN32*/

int cmumps_io_prop_err_info(int ierr);

int cmumps_io_build_err_str(int errnum, int cmumps_err,const char* desc,char* buf,int size);

int MUMPS_CALL cmumps_low_level_init_err_str(char* err_str,int* dim);

