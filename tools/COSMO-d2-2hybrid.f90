program COSMO2hybrid
  use eccodes
  implicit none
  integer            :: ifile,igrib,iret, ofile, ogrib
  character(len=2048) :: val, ifilename, ofilename
  integer            :: grib_count

  real, dimension(132) :: pv

  real, dimension(66), parameter :: &
      &a_half= (/0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,&
      &0.000000, 0.000000, 0.000000, 0.003927, 0.022247, 0.043102, 0.066132, 0.090987, 0.117420, 0.145223, 0.174177, 0.204110,&
      &0.234850, 0.266263, 0.298182, 0.330452, 0.362956, 0.395523, 0.428035, 0.460329, 0.492252, 0.523716, 0.554593, 0.584754,&
      &0.614120, 0.642555, 0.669962, 0.696367, 0.721848, 0.746336, 0.769695, 0.791773, 0.812485, 0.831776, 0.849645, 0.866104,&
      &0.881204, 0.894992, 0.907538, 0.918910, 0.929188, 0.938445, 0.946768, 0.954230, 0.960915, 0.966892, 0.972234, 0.977006,&
      &0.981282, 0.985123, 0.988602, 0.991777, 0.994710, 0.997446, 1.000000/)

  real, dimension(66), parameter :: &
      &b_half = (/4033.8914, 4701.9294, 5451.2650, 6292.3854, 7227.4741, 8266.4106, 9409.7840, 10666.5930, 12035.8470, 13523.8756,&
      &15128.8780, 16860.0323, 18726.5646, 20749.6305, 22487.5720, 22867.6100, 23070.7117, 23120.9218, 23040.1217, 22842.9806,&
      &22539.6543, 22143.5493, 21662.9437, 21106.9771, 20480.3718, 19792.5008, 19050.7574, 18260.4317, 17431.2006, 16568.9728,&
      &15683.3651, 14783.7569, 13873.2668, 12960.0596, 12052.0322, 11154.0865, 10275.3581, 9422.9939, 8592.0471, 7774.5700,&
      &6975.1707, 6206.3803, 5480.9198, 4807.3109, 4190.2981, 3631.3059, 3129.6300, 2682.8609, 2288.0615, 1941.3810, 1638.6812,&
      &1375.7599, 1148.3797, 952.1640, 783.1738, 637.8655, 513.2906, 406.8789, 316.3900, 239.7650, 175.2500, 121.1948, 76.4448,&
      &40.3389, 13.2902, 0.0000 /)
  
  !! a_half
  pv(1:66) = a_half(1:66)
  !! b_half
  pv(67:132) = b_half(1:66)
 
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
!!    call codes_get(igrib, "modelName", val)
    call codes_get(igrib, "nlev", iret)
    if (iret /= 66 ) then
        write (*,*) "Failed to process message ", grib_count
!        write (*,*) "modelName = '", trim(val), "' nlev = ", iret
        write (*,*) "nlev=66 expected"
        call codes_dump(igrib)
        call exit(123)
    endif


    call codes_get(igrib, "typeOfLevel", val)
    !! FOR some reason crashes without clone with eccodes from UBUNTU 18.04
!     if (val == "generalVerticalLayer") then

    call codes_set(igrib,'typeOfLevel','hybridLayer')
    call codes_clone(igrib,ogrib)

    if (val == "generalVertical") then
        !call codes_set(ogrib,'typeOfLevel','hybridLevel')
        call codes_set(ogrib,'typeOfSecondFixedSurface',255)
        !      typeOfSecondFixedSurface=255;
    elseif (val /= "generalVerticalLayer") then
       write (*,*) "Failed to process message ", grib_count
       write (*,*) "typeOfLevel = '", trim(val), "', can be one of generalVerticalLayer or generalVertical"
       call exit(124)
    endif

    call codes_set(ogrib,'PVPresent',1)
    call codes_set(ogrib,'pv',pv)

    !TLossless packing does no harm
    call codes_set(ogrib,'packingType', 'grid_ccsds')
 
    call codes_write(ogrib,ofile)
    call codes_release(ogrib)
    call codes_release(igrib)
  end do
 
 
  call codes_close_file(ifile)
  call codes_close_file(ofile)
 
end program COSMO2hybrid

!
!      65 1.000000   0.0000 1013.2500
!      64 0.997446  13.2902 1010.7955
!      63 0.994710  40.3389 1008.2932
!      62 0.991777  76.4448 1005.6823
!      61 0.988602 121.1948 1002.9132
!      60 0.985123 175.2500 999.9287
!      59 0.981282 239.7650 996.6821
!      58 0.977006 316.3900 993.1155
!      57 0.972234 406.8789 989.1845
!      56 0.966892 513.2906 984.8360
!      55 0.960915 637.8655 980.0259
!      54 0.954230 783.1738 974.7052
!      53 0.946768 952.1640 968.8339
!      52 0.938445 1148.3797 962.3628
!      51 0.929188 1375.7599 955.2573
!      50 0.918910 1638.6812 947.4728
!      49 0.907538 1941.3810 938.9770
!      48 0.894992 2288.0615 929.7315
!      47 0.881204 2682.8609 919.7091
!      46 0.866104 3129.6300 908.8761
!      45 0.849645 3631.3059 897.2154
!      44 0.831776 4190.2981 884.7001
!      43 0.812485 4807.3109 871.3235
!      42 0.791773 5480.9198 857.0736
!      41 0.769695 6206.3803 841.9572
!      40 0.746336 6975.1707 825.9770
!      39 0.721848 7774.5700 809.1577
!      38 0.696367 8592.0471 791.5142
!      37 0.669962 9422.9939 773.0687
!      36 0.642555 10275.3581 753.8223
!      35 0.614120 11154.0865 733.7984
!      34 0.584754 12052.0322 713.0226
!      33 0.554593 12960.0596 691.5419
!      32 0.523716 13873.2668 669.3880
!      31 0.492252 14783.7569 646.6118
!      30 0.460329 15683.3651 623.2622
!      29 0.428035 16568.9728 599.3965
!      28 0.395523 17431.2006 575.0753
!      27 0.362956 18260.4317 550.3693
!      26 0.330452 19050.7574 525.3378
!      25 0.298182 19792.5008 500.0582
!      24 0.266263 20480.3718 474.5950
!      23 0.234850 21106.9771 449.0317
!      22 0.204110 21662.9437 423.4438
!      21 0.174177 22143.5493 397.9201
!      20 0.145223 22539.6543 372.5435
!      19 0.117420 22842.9806 347.4057
!      18 0.090987 23040.1217 322.5942
!      17 0.066132 23120.9218 298.2179
!      16 0.043102 23070.7117 274.3797
!      15 0.022247 22867.6100 251.2175
!      14 0.003927 22487.5720 228.8544
!      13 0.000000 20749.6305 207.4963
!      12 0.000000 18726.5646 187.2656
!      11 0.000000 16860.0323 168.6003
!      10 0.000000 15128.8780 151.2888
!       9 0.000000 13523.8756 135.2388
!       8 0.000000 12035.8470 120.3585
!       7 0.000000 10666.5930 106.6659
!       6 0.000000 9409.7840  94.0978
!       5 0.000000 8266.4106  82.6641
!       4 0.000000 7227.4741  72.2747
!       3 0.000000 6292.3854  62.9239
!       2 0.000000 5451.2650  54.5126
!       1 0.000000 4701.9294  47.0193
!       0 0.000000 4033.8914  40.3389
!
!
!
!
