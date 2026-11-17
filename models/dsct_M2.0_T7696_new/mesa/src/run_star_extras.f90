! Extra profile columns for the nonlinear mode-coupling pipeline, and a
! stop condition on surface gravity.

module run_star_extras

   use star_lib
   use star_def
   use const_def
   use math_lib
   use utils_lib, only: mesa_error
   use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_gamma1, i_gamma3

   implicit none

   ! log g fixes R, hence E_star = GM^2/R and the frequency scale
   ! sqrt(GM/R^3); Teff enters nothing downstream. So the target is log g.
   real(dp), parameter :: default_target_logg = 3.900d0
   real(dp), parameter :: dense_Teff = 8000d0    ! K; below this, save a profile every step
   real(dp), parameter :: zams_h1_drop = 0.02d0  ! central H drop marking the end of the pre-MS

   logical :: warned_gamma3 = .false.

contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      s% extras_startup => extras_startup
      s% extras_start_step => extras_start_step
      s% extras_check_model => extras_check_model
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns

      s% how_many_extra_history_header_items => how_many_extra_history_header_items
      s% data_for_extra_history_header_items => data_for_extra_history_header_items
      s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
      s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

   end subroutine extras_controls


   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      warned_gamma3 = .false.
      write(*,'(a,f9.4)') 'dsct_M2.0: will stop at first post-ZAMS model with log g <= ', target_logg(s)
   end subroutine extras_startup


   real(dp) function target_logg(s)
      type (star_info), pointer :: s
      if (s% x_ctrl(1) > 0d0) then
         target_logg = s% x_ctrl(1)
      else
         target_logg = default_target_logg
      end if
   end function target_logg


   ! Photospheric, matching what select_model.py recomputes from the header.
   real(dp) function photosphere_logg(s)
      type (star_info), pointer :: s
      photosphere_logg = safe_log10(standard_cgrav * s% m(1) / (s% photosphere_r * Rsun)**2)
   end function photosphere_logg


   integer function extras_start_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_start_step = 0

      ! Dense sampling near the target so the final choice is not limited to
      ! the terminating step.
      if (is_post_zams(s) .and. s% Teff < dense_Teff) then
         s% profile_interval = 1
         s% max_years_for_timestep = 5d5
      end if
   end function extras_start_step


   ! Distinguishes the main sequence from the pre-MS, which is also cool and
   ! low-gravity and would otherwise satisfy the stop condition immediately.
   logical function is_post_zams(s)
      type (star_info), pointer :: s
      is_post_zams = (s% center_h1 < 1d0 - s% initial_y - s% initial_z - zams_h1_drop)
   end function is_post_zams


   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going

      if (is_post_zams(s) .and. photosphere_logg(s) <= target_logg(s)) then
         write(*,'(a)') ''
         write(*,'(a)') 'reached the target model'
         write(*,'(a,f9.4)')  '  log g      = ', photosphere_logg(s)
         write(*,'(a,f9.1)')  '  Teff       = ', s% Teff
         write(*,'(a,f9.4)')  '  R / Rsun   = ', s% photosphere_r
         write(*,'(a,f9.4)')  '  log L/Lsun = ', safe_log10(s% photosphere_L)
         write(*,'(a,f9.4)')  '  center h1  = ', s% center_h1
         write(*,'(a)') ''
         s% need_to_save_profiles_now = .true.
         s% save_profiles_model_priority = 10
         extras_check_model = terminate
         termination_code_str(t_xtra1) = 'reached target log g'
         s% termination_code = t_xtra1
         return
      end if

      if (extras_check_model == terminate) s% termination_code = t_extras_check_model
   end function extras_check_model


   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 2
   end function how_many_extra_history_columns


   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: R_cm
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      R_cm = s% photosphere_r * Rsun

      names(1) = 'dyn_freq_cpd'
      vals(1) = sqrt(standard_cgrav * s% m(1) / R_cm**3) * 86400d0 / (2d0*pi)

      names(2) = 'detuning_cut_cpd'
      vals(2) = 0.15d0 * vals(1)
   end subroutine data_for_extra_history_columns


   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 5
   end function how_many_extra_profile_columns


   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz,n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      real(dp) :: res(num_eos_basic_results)
      real(dp) :: d_dlnd(num_eos_basic_results)
      real(dp) :: d_dlnT(num_eos_basic_results)
      real(dp), allocatable :: d_dxa(:,:)
      real(dp) :: Gamma3m1, Gamma3m1_eos, L_floor, t_cum, logRho, logT

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (n /= 5) call mesa_error(__FILE__, __LINE__, 'data_for_extra_profile_columns: expected 5')

      allocate(d_dxa(num_eos_d_dxa_results, s% species))

      names(1) = 'dGamma1_dlnRho_s'   ! at constant entropy
      names(2) = 'dGamma1_dlnRho_T'
      names(3) = 'dGamma1_dlnT_Rho'
      names(4) = 't_thermal'          ! s, thermal time above radius r
      names(5) = 'omega_conv'         ! rad/s

      ! Guards the thermal-timescale integral where L passes through zero.
      L_floor = 1d-6 * max(s% L(1), 1d0)
      t_cum = 0d0

      do k = 1, nz

         ! star_get_eos is the wrapper MESA itself uses, so this reproduces
         ! the exact eos blend that built the model.
         logRho = s% lnd(k) / ln10
         logT = s% lnT(k) / ln10
         call star_get_eos(id, k, s% xa(:,k), &
              s% rho(k), logRho, s% T(k), logT, &
              res, d_dlnd, d_dlnT, d_dxa, ierr)
         if (ierr /= 0) then
            deallocate(d_dxa)
            return
         end if

         ! Gamma_3 - 1 = dlnT/dlnRho|_s, taken as grada*Gamma_1 because
         ! res(i_gamma3) holds Gamma_3 itself, contrary to the eos_def comment.
         Gamma3m1 = s% grada(k) * s% gamma1(k)
         Gamma3m1_eos = res(i_gamma3) - 1d0
         if (.not. warned_gamma3 .and. k == nz/2) then
            if (abs(Gamma3m1 - Gamma3m1_eos) > 1d-6 * max(abs(Gamma3m1), 1d-99)) then
               write(*,*) 'WARNING: Gamma_3-1 mismatch,', Gamma3m1, 'vs', Gamma3m1_eos
            end if
            warned_gamma3 = .true.
         end if

         ! dG1/dlnRho|_s = dG1/dlnRho|_T + (Gamma_3 - 1) * dG1/dlnT|_rho
         vals(k,2) = d_dlnd(i_gamma1)
         vals(k,3) = d_dlnT(i_gamma1)
         vals(k,1) = vals(k,2) + Gamma3m1 * vals(k,3)

         ! int_r^R (T c_P / L) dm; k = 1 is the surface, so this runs inward.
         t_cum = t_cum + s% T(k) * s% Cp(k) * s% dm(k) / max(s% L(k), L_floor)
         vals(k,4) = t_cum

         if (s% conv_vel(k) > 0d0 .and. s% scale_height(k) > 0d0) then
            vals(k,5) = s% conv_vel(k) / (s% mixing_length_alpha * s% scale_height(k))
         else
            vals(k,5) = 0d0
         end if

      end do

      deallocate(d_dxa)
   end subroutine data_for_extra_profile_columns


   integer function how_many_extra_history_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_header_items = 0
   end function how_many_extra_history_header_items


   subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return
   end subroutine data_for_extra_history_header_items


   integer function how_many_extra_profile_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_header_items = 3
   end function how_many_extra_profile_header_items


   subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      real(dp) :: R_cm
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      R_cm = s% photosphere_r * Rsun
      names(1) = 'mixing_length_alpha'
      vals(1) = s% mixing_length_alpha
      names(2) = 'E_star_erg'                         ! G M^2 / R
      vals(2) = standard_cgrav * s% m(1)**2 / R_cm
      names(3) = 'dyn_freq_rad_per_s'                 ! sqrt(G M / R^3)
      vals(3) = sqrt(standard_cgrav * s% m(1) / R_cm**3)
   end subroutine data_for_extra_profile_header_items


   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going
      if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
   end function extras_finish_step


   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      write(*,'(a)') 'run select_model.py to pick the profile and its .GYRE file'
   end subroutine extras_after_evolve

end module run_star_extras
