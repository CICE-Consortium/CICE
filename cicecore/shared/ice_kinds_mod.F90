!=======================================================================

! Defines variable precision for all common data types
! Code originally based on kinds_mod.F in POP
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
! 2006: ECH converted to free source form (F90)

      module ice_kinds_mod

!=======================================================================

      use icepack_intfc, only: char_len  => icepack_char_len
      use icepack_intfc, only: char_len_long => icepack_char_len_long
      use icepack_intfc, only: log_kind  => icepack_log_kind
      use icepack_intfc, only: int_kind  => icepack_int_kind
      use icepack_intfc, only: real_kind => icepack_real_kind
      use icepack_intfc, only: dbl_kind  => icepack_dbl_kind
      use icepack_intfc, only: r16_kind  => icepack_r16_kind

      implicit none
      public

!=======================================================================

      end module ice_kinds_mod

!=======================================================================
