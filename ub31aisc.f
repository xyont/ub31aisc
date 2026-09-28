!     ==================================================================
!     codecheck_aisc.f
!     AISC 360-16 (LRFD) Structural Steel Code Checking
!     ==================================================================

!     ------------------------------------------------------------------
!     1. calc_aisc_properties
!        Menghitung properti penampang + klasifikasi Compact/Non-Compact/Slender
!     ------------------------------------------------------------------
      subroutine calc_aisc_properties(sect_type, dims, E, nu, fy,
     &  A, Iyy, Izz, Sx, Sy, Zx, Zy, rts, ho, J, Cw, &
     &  rx, ry, lambda_f, lambda_w, lambda_pf, lambda_rf, &
     &  lambda_pw, lambda_rw, is_compact, is_slender, is_compact_web)
      implicit none

      real*8,  intent(in)  :: sect_type, dims(6), E, nu, fy
      real*8,  intent(out) :: A, Iyy, Izz, Sx, Sy, Zx, Zy
      real*8,  intent(out) :: rts, ho, J, Cw
      real*8,  intent(out) :: rx, ry
      real*8,  intent(out) :: lambda_f, lambda_w
      real*8,  intent(out) :: lambda_pf, lambda_rf
      real*8,  intent(out) :: lambda_pw, lambda_rw
      integer, intent(out) :: is_compact, is_slender, is_compact_web

      real*8  :: b, h, tf, tw, bt, tt, bb, tb, ro, t, ri
      real*8  :: pi, eps_aisc, c_f, c_w, cf_tf, cw_tw
      integer :: st

      pi = 4.d0 * datan(1.d0)
      st = nint(sect_type)

!     --- Parameter AISC ---
!     E  = modulus elastisitas baja (200e9 Pa untuk baja struktural)
!     fy = tegangan leleh (mis. 250e6, 355e6 Pa)
!
!     Batas kelangsingan AISC Table B4.1a (sayap profil I, tekan):
!       lambda_pf = 0.38 * sqrt(E/fy)
!       lambda_rf = 1.00 * sqrt(E/fy)
!     Batas kelangsingan AISC Table B4.1b (web profil I, lentur):
!       lambda_pw = 3.76 * sqrt(E/fy)      (compact, plastis)
!       lambda_rw = 5.70 * sqrt(E/fy)      (non-compact, inelastis)
!
!     Untuk profil non-I, konstanta berbeda -- lihat Table B4.1a/B4.1b.

      eps_aisc = dsqrt(E / max(fy, 1.0d6))

!     Default nilai-nilai (hindari uninitialized)
      A = 0.01d0; Iyy = 1.0d-5; Izz = 1.0d-5
      Sx = 1.0d-4; Sy = 1.0d-4
      Zx = 1.0d-4; Zy = 1.0d-4
      rts = 0.1d0; ho = 0.2d0
      J = 1.0d-6;  Cw = 0.0d0
      is_compact = 1; is_slender = 0; is_compact_web = 1

!     >>> Di sini tempat menyalin perhitungan geometri A, Iyy, Izz, S, Z
!     >>> dari calc_ec3_properties (bagian geometri dapat dipakai ulang!)
!     >>> Cukup tambahkan: rts, ho, J, Cw, dan batas kelangsingan AISC.

      if (st .eq. 3) then
!        --- I-Beam / W-Shape ---
!        dims: h, b_top, t_f_top, b_bot, t_f_bot, t_w
         h   = dims(1); b_t = dims(2); t_t = dims(3)
         b_b = dims(4); t_b = dims(5); t_w = dims(6)
         b  = 0.5d0*(b_t + b_b)
         tf = 0.5d0*(t_t + t_b)

!        Properti geometri (sama seperti EC3):
!        A, Iyy, Izz, Sx, Sy, Zx, Zy   <-- isi dari perhitungan EC3

!        Radius girasi efektif untuk LTB (AISC F2-7):
         rts = dsqrt( dsqrt(Iyy*Cw) / max(Sx, 1.0d-12) )
!        Jarak antar titik berat sayap (AISC F2-8):
         ho  = h - tf
!        J dan Cw (sama seperti EC3):
!        J  = (1/3)*(2*b*tf**3 + (h-2*tf)*tw**3)
!        Cw = Izz*(h-tf)**2 / 4

!        --- Klasifikasi AISC (profil I) ---
         lambda_f = (b/2.d0) / max(tf, 1.0d-6)   ! b/2tf sayap
         lambda_w = (h - 2.d0*tf) / max(tw, 1.0d-6)  ! h/tw web

         lambda_pf = 0.38d0 * eps_aisc
         lambda_rf = 1.00d0 * eps_aisc
         lambda_pw = 3.76d0 * eps_aisc
         lambda_rw = 5.70d0 * eps_aisc

         if (lambda_f .le. lambda_pf .and.
     &       lambda_w .le. lambda_pw) then
            is_compact = 1
         elseif (lambda_f .le. lambda_rf .and.
     &           lambda_w .le. lambda_rw) then
            is_compact = 0      ! non-compact
         else
            is_slender = 1
         endif

         if (lambda_w .le. lambda_pw) then
            is_compact_web = 1
         else
            is_compact_web = 0
         endif

      elseif (st .eq. 7) then
!        --- Box / RHS (AISC Table B4.1a: b/t = b/t, h/t) ---
         lambda_f = dims(2) / max(dims(3), 1.0d-6)
         lambda_w = dims(1) / max(dims(4), 1.0d-6)
         lambda_pf = 1.12d0 * eps_aisc
         lambda_rf = 1.40d0 * eps_aisc
         lambda_pw = 2.42d0 * eps_aisc
         lambda_rw = 5.70d0 * eps_aisc
         is_compact = merge(1, 0, lambda_f .le. lambda_pf)
         is_slender = merge(1, 0, lambda_f .gt. lambda_rf)

      else
!        --- Fallback: profil lain ---
         lambda_f   = 0.0d0; lambda_w = 0.0d0
         lambda_pf  = 0.38d0*eps_aisc; lambda_rf = 1.0d0*eps_aisc
         lambda_pw  = 3.76d0*eps_aisc; lambda_rw = 5.70d0*eps_aisc
         is_compact = 0; is_slender = 0; is_compact_web = 0
      endif

      rx = dsqrt(max(1.0d-12, Iyy / max(A, 1.0d-12)))
      ry = dsqrt(max(1.0d-12, Izz / max(A, 1.0d-12)))

      end subroutine calc_aisc_properties


!     ------------------------------------------------------------------
!     2. calc_aisc_member_capacities
!        Menghitung kapasitas LRFD: Pn, Vn, Mn, dan kuat tekuk Pn (kolom)
!     ------------------------------------------------------------------
      subroutine calc_aisc_member_capacities(
     &  L_elem, Lb_in, Kx_in, Ky_in, Cb_in,
     &  fy_in, E_in, nu_in, phi_c_in, phi_b_in, phi_v_in,
     &  sect_type, dims,
     &  Py, Pn_yield, Pn_buck_x, Pn_buck_y, Pn_comp, Pn_tens,
     &  Vn_y, Vn_z, Mn_x, Mn_y, Mn_ltb,
     &  Fcr_x, Fcr_y, Lcx, Lcy, lambda_c_x, lambda_c_y,
     &  is_compact, is_slender)
      implicit none

      real*8,  intent(in)  :: L_elem, Lb_in, Kx_in, Ky_in, Cb_in
      real*8,  intent(in)  :: fy_in, E_in, nu_in
      real*8,  intent(in)  :: phi_c_in, phi_b_in, phi_v_in
      real*8,  intent(in)  :: sect_type, dims(6)
      real*8,  intent(out) :: Py, Pn_yield, Pn_buck_x, Pn_buck_y
      real*8,  intent(out) :: Pn_comp, Pn_tens
      real*8,  intent(out) :: Vn_y, Vn_z, Mn_x, Mn_y, Mn_ltb
      real*8,  intent(out) :: Fcr_x, Fcr_y, Lcx, Lcy
      real*8,  intent(out) :: lambda_c_x, lambda_c_y
      integer, intent(out) :: is_compact, is_slender

      real*8 :: A, Iyy, Izz, Sx, Sy, Zx, Zy, rts, ho, J, Cw
      real*8 :: rx, ry, lambda_f, lambda_w
      real*8 :: lambda_pf, lambda_rf, lambda_pw, lambda_rw
      integer :: is_compact_web

      real*8 :: pi, Lb, Kx, Ky, Cb, fy, E, G
      real*8 :: Lp, Lr, BF, Mp, Mr, Mn_ltb_val
      real*8 :: Fcr, lambda_c, Fe

!     --- Panggil properti ---
      call calc_aisc_properties(sect_type, dims, E_in, nu_in, fy_in,
     &  A, Iyy, Izz, Sx, Sy, Zx, Zy, rts, ho, J, Cw,
     &  rx, ry, lambda_f, lambda_w, lambda_pf, lambda_rf,
     &  lambda_pw, lambda_rw, is_compact, is_slender, is_compact_web)

      pi = 4.d0 * datan(1.d0)
      Lb = max(0.01d0, Lb_in)
      Kx = max(0.1d0, Kx_in)
      Ky = max(0.1d0, Ky_in)
      Cb = max(1.0d0, Cb_in)
      fy = max(1.0d6, fy_in)
      E  = max(1.0d9, E_in)
      G  = E / (2.d0*(1.d0 + nu_in))

!     --- Kuat Leleh ---
      Py      = fy * A
      Pn_tens = fy * A              ! AISC D2-1
      Pn_yield= fy * A              ! untuk komponen tekan leleh

!     --- Kuat Tekan Kolom (AISC E3) ---
      Lcx = Kx * L_elem
      Lcy = Ky * L_elem
      lambda_c_x = (Lcx / max(rx,1.d-6)) / pi * dsqrt(fy/E)
      lambda_c_y = (Lcy / max(ry,1.d-6)) / pi * dsqrt(fy/E)

      if (lambda_c_x .le. 1.5d0) then
         Fcr_x = (0.658d0**(lambda_c_x**2)) * fy
      else
         Fcr_x = (0.877d0 / (lambda_c_x**2)) * fy
      endif
      if (lambda_c_y .le. 1.5d0) then
         Fcr_y = (0.658d0**(lambda_c_y**2)) * fy
      else
         Fcr_y = (0.877d0 / (lambda_c_y**2)) * fy
      endif

      Pn_buck_x = Fcr_x * A
      Pn_buck_y = Fcr_y * A
      Pn_comp   = min(Pn_buck_x, Pn_buck_y)   ! kritis terhadap sumbu lemah

!     --- Kuat Geser (AISC G2) ---
      Vn_y = 0.6d0 * fy * (2.d0 * Sx)         ! aproksimasi Aw_y
      Vn_z = 0.6d0 * fy * (2.d0 * Sy)         ! aproksimasi Aw_z
!     Untuk presisi, gunakan Aw = tinggi web × tebal web

!     --- Kuat Lentur (AISC F2) ---
      Mp = fy * Zx                            ! momen plastis
      Mr = 0.7d0 * fy * Sx                    ! momen batas inelastis

!     Panjang tanpa bracing pembatas (AISC F2-5 & F2-6):
      Lp = 1.76d0 * ry * dsqrt(E / fy)
      Lr = 1.95d0 * rts * (E / (0.7d0*fy)) *
     &     dsqrt( J/(Sx*ho) +
     &       dsqrt( (J/(Sx*ho))**2 +
     &         6.76d0 * (0.7d0*fy/E)**2 ) )

      if (is_compact .eq. 1) then
         if (Lb .le. Lp) then
!           Zona plastis penuh
            Mn_x = Mp
         elseif (Lb .le. Lr) then
!           Zona inelastis LTB (AISC F2-2)
            BF = (Mp - Mr) / max(Lr - Lp, 1.d-9)
            Mn_x = min( Mp, Cb*(Mp - BF*(Lb - Lp)) )
         else
!           Zona elastis LTB (AISC F2-3 & F2-4)
            Fe = (Cb * pi**2 * E / (Lb/rts)**2) *
     &           dsqrt( 1.d0 + 0.078d0 *
     &                  (J*Cw/(Sx*ho))**2 / (Lb/rts)**2 )
!           Koreksi: gunakan rumus F2-4 untuk Fe
            Fe = (Cb * pi**2 * E) / (Lb/rts)**2 *
     &           dsqrt( 1.d0 + 0.078d0 *
     &                  (J / (Sx*ho)) * (Lb/rts)**2 )
            Mn_x = min( Mp, 0.9d0 * Fe * Sx )
         endif
      else
!        Penampang non-compact: interpolasi Mp - Mr
         Mn_x = Mr + (Mp - Mr) *
     &          (lambda_rf - lambda_f) /
     &          max(lambda_rf - lambda_pf, 1.d-9)
         Mn_x = min(Mn_x, Mp)
      endif
      Mn_ltb = Mn_x
      Mn_y   = fy * Zy                        ! lentur sumbu lemah

      end subroutine calc_aisc_member_capacities


!     ------------------------------------------------------------------
!     3. calc_aisc_station_check
!        Cek rasio pemanfaatan (UC) per stasiun, format sama seperti EC3
!     ------------------------------------------------------------------
      subroutine calc_aisc_station_check(fx, vy, vz, tx, my, mz,
     &  Py, Pn_tens, Pn_comp, Vn_y, Vn_z, Mn_x, Mn_y, Mn_ltb,
     &  phi_c, phi_b, phi_v, is_compact,
     &  uc_tot, uc_ax, uc_sh, uc_bnd, uc_buck, uc_ltb, uc_stab,
     &  gov_mode)
      implicit none

      real*8,  intent(in)  :: fx, vy, vz, tx, my, mz
      real*8,  intent(in)  :: Py, Pn_tens, Pn_comp
      real*8,  intent(in)  :: Vn_y, Vn_z, Mn_x, Mn_y, Mn_ltb
      real*8,  intent(in)  :: phi_c, phi_b, phi_v
      integer, intent(in)  :: is_compact
      real*8,  intent(out) :: uc_tot, uc_ax, uc_sh, uc_bnd
      real*8,  intent(out) :: uc_buck, uc_ltb, uc_stab
      character*20, intent(out) :: gov_mode

      real*8 :: Pr, Mr_x, Mr_y, phiPn_c, phiPn_t, phiVn_y, phiVn_z
      real*8 :: phiMn_x, phiMn_y, phiMn_ltb
      real*8 :: uc_H1a, uc_H1b, ratio_pr_pc
      real*8 :: uc_bnd_y, uc_bnd_z, uc_ax_t, uc_ax_c
      real*8 :: phi_Pn_buck

!     --- Kapasitas terfaktor (LRFD) ---
      phiPn_t  = phi_c * Pn_tens
      phiPn_c  = phi_c * Pn_comp
      phiVn_y  = phi_v * Vn_y
      phiVn_z  = phi_v * Vn_z
      phiMn_x  = phi_b * Mn_x
      phiMn_y  = phi_b * Mn_y
      phiMn_ltb= phi_b * Mn_ltb

!     --- 1. Aksial ---
      uc_ax_t = 0.d0; uc_ax_c = 0.d0
      if (fx .gt. 0.d0) then
         uc_ax_t = fx / max(phiPn_t, 1.d-6)
         uc_ax   = uc_ax_t
      else
         uc_ax_c = dabs(fx) / max(phiPn_c, 1.d-6)
         uc_ax   = uc_ax_c
      endif

!     --- 2. Geser ---
      uc_sh = max( dabs(vy)/max(phiVn_y,1.d-6),
     &             dabs(vz)/max(phiVn_z,1.d-6) )

!     --- 3. Lentur ---
      uc_bnd_y = dabs(my) / max(phiMn_x, 1.d-6)
      uc_bnd_z = dabs(mz) / max(phiMn_y, 1.d-6)
      uc_bnd   = uc_bnd_y + uc_bnd_z

!     --- 4. Tekuk Kolom (Compression buckling) ---
      if (fx .lt. 0.d0) then
         uc_buck = dabs(fx) / max(phiPn_c, 1.d-6)
      else
         uc_buck = 0.d0
      endif

!     --- 5. LTB ---
      uc_ltb = dabs(my) / max(phiMn_ltb, 1.d-6)

!     --- 6. Interaksi Tekan + Lentur (AISC H1-1a / H1-1b) ---
      if (fx .lt. 0.d0) then
         Pr      = dabs(fx)
         phi_Pn_buck = phiPn_c
         ratio_pr_pc = Pr / max(phi_Pn_buck, 1.d-6)

         Mr_x = dabs(my) / max(phiMn_ltb, 1.d-6)
         Mr_y = dabs(mz) / max(phiMn_y, 1.d-6)

         if (ratio_pr_pc .ge. 0.2d0) then
!           AISC Eq. H1-1a
            uc_H1a = ratio_pr_pc +
     &               (8.d0/9.d0)*(Mr_x + Mr_y)
            uc_stab = uc_H1a
         else
!           AISC Eq. H1-1b
            uc_H1b = 0.5d0*ratio_pr_pc + (Mr_x + Mr_y)
            uc_stab = uc_H1b
         endif
      else
         uc_stab = uc_ax_t + uc_bnd
      endif

!     --- 7. Mode penentu ---
      uc_tot = max(uc_ax, uc_sh, uc_bnd, uc_buck, uc_ltb, uc_stab)

      if (uc_tot .eq. uc_stab) then
         if (fx .lt. 0.d0) then
            gov_mode = 'AISC-H1-1'
         else
            gov_mode = 'AISC-Tens+Bend'
         endif
      elseif (uc_tot .eq. uc_ltb) then
         gov_mode = 'AISC-LTB'
      elseif (uc_tot .eq. uc_buck) then
         gov_mode = 'AISC-Compression'
      elseif (uc_tot .eq. uc_bnd) then
         gov_mode = 'AISC-Bending'
      elseif (uc_tot .eq. uc_sh) then
         gov_mode = 'AISC-Shear'
      elseif (uc_tot .eq. uc_ax) then
         if (fx .gt. 0.d0) then
            gov_mode = 'AISC-Tension'
         else
            gov_mode = 'AISC-Compression'
         endif
      else
         gov_mode = 'AISC-General'
      endif

      end subroutine calc_aisc_station_check


!     ------------------------------------------------------------------
!     4. calc_aisc_member_check
!        Wrapper utama, format paralel dengan calc_ec3_member_check
!     ------------------------------------------------------------------
      subroutine calc_aisc_member_check(fx, vy, vz, tx, my, mz,
     &  L_elem, Lb_in, Kx_in, Ky_in, Cb_in,
     &  fy_in, E_in, nu_in, phi_c_in, phi_b_in, phi_v_in,
     &  sect_type, dims,
     &  uc_tot, uc_ax, uc_sh, uc_bnd, uc_buck, uc_ltb, uc_stab,
     &  gov_mode)
      implicit none

      real*8,  intent(in)  :: fx, vy, vz, tx, my, mz
      real*8,  intent(in)  :: L_elem, Lb_in, Kx_in, Ky_in, Cb_in
      real*8,  intent(in)  :: fy_in, E_in, nu_in
      real*8,  intent(in)  :: phi_c_in, phi_b_in, phi_v_in
      real*8,  intent(in)  :: sect_type, dims(6)
      real*8,  intent(out) :: uc_tot, uc_ax, uc_sh, uc_bnd
      real*8,  intent(out) :: uc_buck, uc_ltb, uc_stab
      character*20, intent(out) :: gov_mode

      real*8  :: Py, Pn_tens, Pn_comp
      real*8  :: Vn_y, Vn_z, Mn_x, Mn_y, Mn_ltb
      real*8  :: Fcr_x, Fcr_y, Lcx, Lcy, lambda_c_x, lambda_c_y
      real*8  :: phi_c, phi_b, phi_v
      integer :: is_compact, is_slender

!     Set faktor resistensi default AISC LRFD:
      phi_c = 0.90d0
      phi_b = 0.90d0
      phi_v = 1.00d0

      call calc_aisc_member_capacities(L_elem, Lb_in, Kx_in, Ky_in,
     &  Cb_in, fy_in, E_in, nu_in, phi_c, phi_b, phi_v,
     &  sect_type, dims,
     &  Py, Pn_tens, Pn_comp, Pn_tens, Pn_comp, Pn_tens,
     &  Vn_y, Vn_z, Mn_x, Mn_y, Mn_ltb,
     &  Fcr_x, Fcr_y, Lcx, Lcy, lambda_c_x, lambda_c_y,
     &  is_compact, is_slender)

      call calc_aisc_station_check(fx, vy, vz, tx, my, mz,
     &  Py, Pn_tens, Pn_comp, Vn_y, Vn_z, Mn_x, Mn_y, Mn_ltb,
     &  phi_c, phi_b, phi_v, is_compact,
     &  uc_tot, uc_ax, uc_sh, uc_bnd, uc_buck, uc_ltb, uc_stab,
     &  gov_mode)

      end subroutine calc_aisc_member_check