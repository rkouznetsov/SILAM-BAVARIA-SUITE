program ICON2hybrid
  use eccodes
  implicit none
  integer            :: ifile,igrib,iret, ofile, ogrib, iTmp, jTmp
  character(len=2048) :: val, ifilename, ofilename
  integer            :: grib_count

  real, dimension(132) :: pv
  real, dimension(65), parameter :: &
!      &a_full=(/ 3863.6357, 4468.8536, 5123.1364, 5823.0421, 6565.7566, 7347.9516, 8155.2128, 8983.9827, 9841.9193, 10708.5316,&
!      &11552.2356, 12445.7914, 13354.6124, 14290.3118, 15251.6678, 16243.4213, 17287.7951, 18387.6791, 19538.9653, 20738.8273, &
!      &21979.8887, 23250.5824, 24538.4088, 25834.3468, 27135.9166, 28442.7547, 29751.3182, 31053.5988, 32330.6430, 33555.1790, &
!      &34702.8267, 35750.6465, 36671.5183, 37429.4103, 37987.0413, 38306.9681, 38357.0492, 38117.4374, 37579.5862, 36744.4370, &
!      &35617.1172, 34208.3297, 32534.5331, 30620.8910, 28492.9989, 26179.2589, 23720.3665, 21160.0048, 18541.9059, 15922.3695, &
!      &13358.5456, 10900.8107, 8602.9857, 6511.0676, 4665.3171, 3101.7344, 1847.5712, 920.0198, 239.9878, 0.0000/)

        & a_full= (/4693.0636, 6452.2761, 7867.3732, 9258.9305, 10666.2949, 12100.8869, &
        & 13565.6432, 15064.4409, 16599.2110, 18172.8825, 19791.4311, 21459.1711, &
        & 23174.1517, 24930.1085, 26721.1341, 28538.8740, 30373.9246, 32220.0970, &
        & 34072.9567, 35814.7505, 36877.1127, 37812.2405, 38619.4158, 39297.5178, &
        & 39845.9003, 40262.7620, 40546.6108, 40697.2564, 40713.5083, 40595.2197, &
        & 40344.6280, 39964.2068, 39456.2582, 38823.1210, 38067.4073, 37192.4698, &
        & 36204.6360, 35111.7073, 33920.2627, 32637.4514, 31271.6695, 29830.4828, &
        & 28321.0715, 26752.8829, 25137.4671, 23487.4315, 21815.2142, 20131.9224, &
        & 18448.3846, 16775.9158, 15124.9610, 13505.5719, 11928.2603, 10404.0587, &
        & 8944.2910, 7559.3756, 6258.7278, 5051.7235, 3947.5999, 2955.5123, 2084.5904, &
        & 1344.5310, 740.7357, 232.5469, 0.0000 /)
     
  real, dimension(65), parameter :: &
 !     &b_full =(/0.001117, 0.000957, 0.000763, 0.000550, 0.000329, 0.000109, 0.000000, 0.000000, 0.000000, 0.000145, 0.000714,&
 !     &0.000935, 0.001116, 0.001245, 0.001457, 0.001850, 0.002361, 0.002967, 0.003728, 0.004699, 0.005987, 0.007734, 0.010085,&
 !     &0.013144, 0.016942, 0.021486, 0.026824, 0.033059, 0.040408, 0.049170, 0.059613, 0.071989, 0.086596, 0.103824, 0.124076,&
 !     &0.147344, 0.173436, 0.202420, 0.234315, 0.269034, 0.306451, 0.346377, 0.388572, 0.432698, 0.478413, 0.525349, 0.573009,&
 !     &0.620866, 0.668394, 0.714931, 0.759808, 0.802429, 0.842155, 0.878428, 0.910738, 0.938607, 0.961627, 0.979459, 0.992647,&
 !     &0.998825/)
      
        & b_full= (/0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
        & 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  &
        & 0.000000, 0.000000, 0.000000, 0.000000, 0.001162, 0.009235, 0.018579, 0.029176,  &
        & 0.041019, 0.054100, 0.068426, 0.084003, 0.100819, 0.118872, 0.138151, 0.158621,  &
        & 0.180245, 0.202989, 0.226816, 0.251686, 0.277554, 0.304345, 0.331970, 0.360351,  &
        & 0.389408, 0.419048, 0.449188, 0.479752, 0.510635, 0.541705, 0.572816, 0.603826,  &
        & 0.634607, 0.665032, 0.694973, 0.724309, 0.752923, 0.780694, 0.807491, 0.833178,  &
        & 0.857624, 0.880706, 0.902296, 0.922267, 0.940481, 0.956792, 0.971029, 0.983031,  &
        & 0.992999, 0.998731/)

  
  !! a_half
  pv(1) = 0
  pv(2:65) = 0.5*(a_full(1:64)+a_full(2:65))
  pv(66) = 0
  !! b_half
  pv(67) = 0
  pv(68:131) = 0.5*(b_full(1:64)+b_full(2:65))
  pv(132) = 1.0
 
  !First, make sure the right number of inputs have been provided
  IF(COMMAND_ARGUMENT_COUNT().NE.2)THEN
    WRITE(*,*)'ERROR, TWO COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
    WRITE(*,*)'Usage: CON2hybrid infile.grib outfile.grib'
    STOP
  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,ifilename)   !first, read in the two values
  CALL GET_COMMAND_ARGUMENT(2,ofilename)

   
  


  call codes_open_file(ifile, trim(ifilename), 'r')
 
  call codes_open_file(ofile, trim(ofilename), 'w')

  ! Loop on all the messages in a file.
  grib_count=0
  do while (.True.)
    call codes_grib_new_from_file(ifile,igrib, iret)
    if (iret == CODES_END_OF_FILE) exit
 
    grib_count=grib_count+1
    !    write(*,*) '-- GRIB N. ',grib_count,' --'

    !! Only for icon levels, other fields pass as_is
    call codes_get(igrib, "typeOfFirstFixedSurface", iret)
    if (iRet /= 150) then
      call codes_write(igrib,ofile)
      call codes_release(igrib)
      cycle
    endif


    !!!! Make sure that it is a right grib...
    call codes_get(igrib, "generatingProcessIdentifier", iTmp)
    call codes_get(igrib, "nlev", jTmp)
    if (jTmp /= 66 .or. 11 /= iTmp) then
        write (*,*) "Failed to process message ", grib_count
        write (*,*) "generatingProcessIdentifier = '", iTmp, "' nlev = ", jTmp
        write (*,*) "modelName=icreu_0.625l60 and nlev=61 expected"
        call codes_dump(igrib)
        call exit(123)
    endif


    call codes_set(igrib,'typeOfLevel','hybridLayer')

    !! FOR some reason crashes without clone with eccodes from UBUNTU 18.04
    call codes_clone(igrib,ogrib)
    call codes_set(ogrib,'PVPresent',1)
    call codes_set(ogrib,'pv',pv)
 
    call codes_write(ogrib,ofile)
    call codes_release(ogrib)
    call codes_release(igrib)
  end do
 
 
  call codes_close_file(ifile)
  call codes_close_file(ofile)
 
end program ICON2hybrid

