module orl

   use scala3, only: nscal,n1,n2,n3

   ! quantita' necessarie per subroutine orlansky ed inflow

   real du_dx1(n2,n3),dv_dx1(n2,n3),dw_dx1(n2,n3)
   real du_dx2(n2,n3),dv_dx2(n2,n3),dw_dx2(n2,n3)

   real du_dy3(n1,n3),dv_dy3(n1,n3),dw_dy3(n1,n3)
   real du_dy4(n1,n3),dv_dy4(n1,n3),dw_dy4(n1,n3)

   real du_dz5(n1,n2),dv_dz5(n1,n2),dw_dz5(n1,n2)
   real du_dz6(n1,n2),dv_dz6(n1,n2),dw_dz6(n1,n2)

   real drho_dx1(nscal,n2,n3),drho_dx2(nscal,n2,n3)
   real drho_dy3(nscal,n1,n3),drho_dy4(nscal,n1,n3)
   real drho_dz5(nscal,n1,n2),drho_dz6(nscal,n1,n2)

   real area_bagnata1,area_bagnata2
   real area_bagnata3,area_bagnata4
   real area_bagnata5,area_bagnata6

   integer index_out1(0:n2+1,0:n3+1),index_out2(0:n2+1,0:n3+1)
   integer index_out3(0:n1+1,0:n3+1),index_out4(0:n1+1,0:n3+1)
   integer index_out5(0:n1+1,0:n2+1),index_out6(0:n1+1,0:n2+1)

   integer index_rho1(0:n2+1,0:n3+1),index_rho2(0:n2+1,0:n3+1)
   integer index_rho3(0:n1+1,0:n3+1),index_rho4(0:n1+1,0:n3+1)
   integer index_rho5(0:n1+1,0:n2+1),index_rho6(0:n1+1,0:n2+1)

   integer infout1,infout2,infout3
   integer infout4,infout5,infout6
      
   integer npn1,npn2,npn3,npn4,npn5,npn6

   common/orlCommon/du_dx1,dv_dx1,dw_dx1, &
      du_dx2,dv_dx2,dw_dx2, &
      du_dy3,dv_dy3,dw_dy3, &
      du_dy4,dv_dy4,dw_dy4, &
      du_dz5,dv_dz5,dw_dz5, &
      du_dz6,dv_dz6,dw_dz6, &
      drho_dx1,drho_dx2, &
      drho_dy3,drho_dy4, &
      drho_dz5,drho_dz6, &
      area_bagnata1,area_bagnata2, &
      area_bagnata3,area_bagnata4, &
      area_bagnata5,area_bagnata6, &
      index_out1,index_out2, &
      index_out3,index_out4, &
      index_out5,index_out6, &
      index_rho1,index_rho2, &
      index_rho3,index_rho4, &
      index_rho5,index_rho6, &
      infout1,infout2,infout3, &
      infout4,infout5,infout6, &
      npn1,npn2,npn3,npn4,npn5,npn6
end module orl

