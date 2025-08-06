module parallel_dynamics

    use debug_flags, only: debug => parallel_dynamics_debug

    implicit none
    private

    public :: init_parallel_dynamics
    public :: advance_parallel_dynamics
    public :: finish_parallel_dynamics
    
    real, dimension (:, :, :, :, :), allocatable :: departure_point_izext, departure_point_iv
    logical, dimension (:, :, :, :, :), allocatable :: departure_point_outside_grid
    
contains
    subroutine init_parallel_dynamics
        implicit none

        ! allocate arrays needed during time advance
        if (debug) write (*, *) 'parallel_dynamics::init_parallel_dynamics::allocate_arrays'
        call allocate_arrays

        ! for every grid point in the 2D (zext, vpa) domain, the velocities in the zext and vpa phase space
        ! are computed in the approximation that the grid points follow their corresponding lowest-order 
        ! characteristic (which has constant particle kinetic energy).

        if (debug) write (*, *) 'parallel_dynamics::init_parallel_dynamics::find_departure_points'
        call find_departure_points

    end subroutine init_parallel_dynamics

    subroutine allocate_arrays
        use parameters_kxky_grids, only: nakx
        use zgrid, only: nzgrid, ntubes
        use vpamu_grids, only: nvpa
        use stella_layouts, only: kymus_lo
        
        implicit none

        if (.not. allocated(departure_point_izext)) allocate(departure_point_izext(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
        if (.not. allocated(departure_point_iv)) allocate(departure_point_iv(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
        if (.not. allocated(departure_point_outside_grid)) allocate(departure_point_outside_grid(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
    end subroutine allocate_arrays

    subroutine find_departure_points
        ! For each point in the (zext, vpa) plane, the departure point is obtained by using the velocity fields to project back
        ! in time. The velocity fields are defined by the following non-dimensionalised expressions: 
        ! vz = vpa * (bhat . grad z)
        ! v_vpa = - mu * (bhat . grad z) dB/dz

        ! Trace back in time using implicit midpoint i.e. dx = dt * U(x-dx/2)
        ! NOTE: EVERYTHING ASSUMES UNIFORM SPACING IN Z AND V GRIDS (seems like the grids are only defined as uniform in the whole code anyway)
        
        use zgrid, only: nzgrid, ntubes, delzed
        use vpamu_grids, only: nvpa, vpa, mu, nmu
        use stella_layouts, only: kymus_lo
        use stella_layouts, only: imu_idx, iky_idx
        use extended_zgrid, only: nsegments, nzed_segment, neigen
        use extended_zgrid, only: map_to_iz_ikx_from_izext
        
        use stella_time, only: code_dt
        use geometry, only: b_dot_grad_z, dbdzed
  
    
        implicit none

        integer :: ikymus, iky, imu
        integer :: ikx, iz, it, iv
        integer :: izext, nz_ext, ie
        integer :: ia = 1  ! flux annulus is not supported
        integer :: ntime
        integer, dimension (:), allocatable :: iz_from_izext, ikx_from_izext
        real, dimension(:), allocatable :: v_zed_from_izext, v_vpa_from_izext
        real :: max_dvzdzed, max_dvvdzed, sub_dt
        
        
        ! v_zed = vpa * (b . grad z)
        ! v_vpa = - mu * (b . grad z) * dBdz
        ! To estimate the upper bound on time-step size, we need to calculate the maximum velocity gradients. Due to the form of
        ! these expressions, we only need to evaluate the gradients in z of the z dependent portion once before beginning the loops
        call find_max_dvdzed(max_dvzdzed, max_dvvdzed)


        ! the characteristic at a given (zext, vpa) value is determined by the particle's energy and mu
        ! energy = v^2 / vths^2 = energy(iv, imu, iz)
        ! the departure point will thus be a function of iv, imu, and izext (i.e., iz, ikx and it)

        ! The indices available are tube number, kx, ky, z, vpa, mu, species. They are grouped in the following manner:
        ! ky, mu and species are independent, and they are indexed by ikymus, which can be reversed. This is looped through
        ! in the outer loop.
        ! kx and z are combined into a larger grid labelled zext. g will be advected along zext and vpa

        do ikymus = kymus_lo%llim_proc, kymus_lo%ulim_proc
            ! Retrieve indices for ky and mu
            imu = imu_idx(kymus_lo, ikymus)
            iky = iky_idx(kymus_lo, ikymus)

            ! Determine required sub time-step and step count
            call find_tstep_count(mu(imu), max_dvzdzed, max_dvvdzed, ntime)
            sub_dt = code_dt/ntime

            print *, "ntime is :" , ntime
            !> FLAG: ok it seems like ntime is just 1 always, so code_dt is the bottleneck...

            do it = 1, ntubes
                do ie = 1, neigen(iky)
                    ! nz_ext is the number of grid points in the extended zed domain
                    nz_ext = nsegments(ie, iky) * nzed_segment + 1
                    ! obtain the mapping to iz and ikx from the extended zed domain
                    allocate(iz_from_izext(nz_ext))
                    allocate(ikx_from_izext(nz_ext))
                    allocate(v_zed_from_izext(nz_ext))
                    allocate(v_vpa_from_izext(nz_ext))

                    call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)
                    
                    ! Note that v_zed is a product of vpa and some function of zed
                    ! v_vpa is only a function of zed
                    ! Therefore it's simpler and less costly to just interpolate in 2 directions separately
                
                    ! Therefore, let's set up the 1D velocity arrays in the zext axis
                    do izext = 1, nz_ext
                        iz = iz_from_izext(izext)
                        v_zed_from_izext(izext) = b_dot_grad_z(ia, iz) ! the factor of vpa will be added later !!!
                        v_vpa_from_izext(izext) = -mu(imu) * b_dot_grad_z(ia, iz) * dbdzed(ia, iz)
                    end do

                    ! Now loop through all grid points (zext, vpa)
                    do iv = 1, nvpa
                        do izext = 1, nz_ext
                            ! Obtain corresponding ikx and iz from initial iv and izext
                            ikx = ikx_from_izext(izext)
                            iz = iz_from_izext(izext)

                            ! Call subroutine to calculate trajectory and find departure point for current iv and izext
                            call calculate_zed_vpa_departure_idx(ntime, sub_dt, nz_ext, izext, iv, &
                                v_zed_from_izext, v_vpa_from_izext, & 
                                departure_point_izext(ikx, iz, it, iv, ikymus), &
                                departure_point_iv(ikx, iz, it, iv, ikymus), &
                                departure_point_outside_grid(ikx, iz, it, iv, ikymus))
                        end do
                    end do
                    
                    deallocate(iz_from_izext, ikx_from_izext, v_zed_from_izext, v_vpa_from_izext)

                end do
            end do

             
        end do


    end subroutine find_departure_points





        ! ##### Assisting subroutines #####!

    subroutine find_max_dvdzed(max_dvzdzed, max_dvvdzed)
        ! The aim is the calculate the maximum gradients in advection velocities to set an upper bound on the sub time-step
        ! This assumes that:
        ! v_zed = vpa * b . grad z 
        ! v_vpa = (b . grad z) dBdz
        ! The factor of mu in v_vpa is dropped and will be included later when evaluating the upper bound of time-step
        ! This is so we don't have to repeat all the calculations below for each mu

        use geometry, only: b_dot_grad_z, dbdzed 
        use zgrid, only: nzgrid, delzed
        use vpamu_grids, only: vpa, nvpa

        real, intent(out) :: max_dvzdzed, max_dvvdzed
        real, dimension(-nzgrid:nzgrid) :: d_vz, d_vv
        integer :: ia = 1 ! No flux annulus


        call get_dzed(nzgrid, delzed, b_dot_grad_z(ia, :), d_vz)
        call get_dzed(nzgrid, delzed, b_dot_grad_z(ia, :) * dbdzed(ia, :), d_vv)

        ! Along the vpa axis, there is only z-dependence of the advection velocity, so simply return the maximum of the z gradient
        max_dvvdzed = maxval(abs(d_vv))

        ! Along the zed axis, the advection velocity is vpa * f(z), so we return the maximum of the gradient in vpa and zed 
        max_dvzdzed = max(maxval(abs(d_vz)) * maxval(abs(vpa)), maxval(abs(b_dot_grad_z(ia, :))))
    end subroutine

    subroutine get_dzed(nz, dz, f, df)
        implicit none

        integer, intent(in) :: nz
        real, dimension(-nz:), intent(in) :: dz, f
        real, dimension(-nz:), intent(out) :: df

        df(-nz + 1:nz - 1) = (f(-nz + 2:) - f(:nz - 2)) / (dz(:nz - 2) + dz(-nz + 1:nz - 1))
        df(-nz) = (f(-nz + 1) - f(nz - 1)) / (dz(-nz) + dz(nz - 1))
        df(nz) = df(-nz)
   	

    end subroutine get_dzed

    subroutine find_tstep_count(mu, max_dvzdzed, max_dvvdzed, number_of_steps)
        ! Here we calculate the number of sub time-steps needed to advance the parallel streaming
        ! Sub time-step is given by the bound dt <= c0/max(grad v), where c0 is an adjustable constant which dictates error 
        ! tolerance, here set to 0.2.

        ! There could be a way to relax c0 to a higher value (e.g. maybe 0.5), although a better way of solving the implicit equation
        ! in each sub-step is needed. Currently it costs more to use a higher c0 because more iterations will be needed for 
        ! convergence (I think...)

        use stella_time, only: code_dt

        real, intent(in) :: mu, max_dvzdzed, max_dvvdzed
        integer, intent(out) :: number_of_steps
        real :: max_dt


        max_dt = 0.1 / max(max_dvzdzed, abs(mu * max_dvvdzed)) 

        number_of_steps = ceiling(code_dt/max_dt)

    end subroutine find_tstep_count

    subroutine calculate_zed_vpa_departure_idx(no_of_steps, dt, nz_ext, initial_izext, initial_ivpa, &
                                v_zed_from_izext, v_vpa_from_izext, & 
                                depart_izext, depart_ivpa, departure_outside_grid)

        use zgrid, only: delzed 
        use vpamu_grids, only: vpa, nvpa, dvpa, vpa_max ! This is is just a real number it seems for dvpa

        integer, intent(in) :: no_of_steps, nz_ext, initial_izext, initial_ivpa
        real, intent(in) :: dt
        real, dimension(:) :: v_zed_from_izext, v_vpa_from_izext
        real, intent(out) :: depart_izext, depart_ivpa
        logical, intent(out) :: departure_outside_grid

        integer :: nt, k
        real :: izext, ivpa, d_izext, d_ivpa, dz
        real, dimension(2) :: dummy_points
        
        dz = delzed(0) ! It is assumed that z-grid is uniform, or this will not work
        departure_outside_grid = .false. ! Starts inside grid

        ! Initialise the transient iz and iv indices, as well as initial guesses for steps in (iz, ivpa)
        izext = initial_izext
        ivpa = initial_ivpa

        d_izext = -dt/2.0/dz * v_zed_from_izext(initial_izext) * vpa(initial_ivpa) 
        d_ivpa = -dt/2.0/dvpa * v_vpa_from_izext(initial_izext) 

        
        timeloop: do nt = 1, no_of_steps
            ! Iterate to solve implicit equation
            ! The final d_izext and d_ivpa is reused as an initial guess, for the next time step.
            do k = 1, 20
                ! Iterate for convergence. Since this is only done once during init, it shouldn't be a computational bottleneck
                
                ! Check if projected departure point has gone outside of grid
                if (int(izext + d_izext) < 1 .or. int(izext + d_izext) >= real(nz_ext)) then
                    departure_outside_grid = .true.
                    ! Particle has exited extended grid, so no need to continue
                    exit timeloop
                end if


                ! v_zed = vpa * v_zed_from_izext(current zext), so we need to interpolate v_zed_from_izext
                ! vpa assumes vpa grid is uniform, and linearly extrapolated beyond the boundaries for current calculation
                dummy_points = [ v_zed_from_izext(int(izext + d_izext)), v_zed_from_izext(int(izext + d_izext)+ 1) ]
                d_izext = -dt/2.0/dz * linear_interpolate(dummy_points, modulo(izext + d_izext, 1.0)) & 
                        * (- vpa_max + dvpa * (ivpa + d_ivpa - 1))
                        

                ! Similarly, v_vpa = v_vpa_from_izext(current zext)
                dummy_points = [ v_vpa_from_izext(int(izext + d_izext)), v_vpa_from_izext(int(izext + d_izext)+ 1) ]
                d_ivpa = -dt/2.0/dvpa * linear_interpolate(dummy_points, modulo(izext + d_izext, 1.0))

            end do

            ! Update izext and ivpa for this time step
            izext = izext + 2.0*d_izext
            ivpa = ivpa + 2.0*d_ivpa

        end do timeloop
        


        ! Return values for indices of departure point
        depart_izext = izext
        depart_ivpa = ivpa

        ! Now check if particle has departed the boundaries of vpa grid
        ! FLAG: need to define departure iv in this case... only relevant for the maxwellian acceleration term
        if (ivpa > nvpa) then
            departure_outside_grid = .true. 
            ivpa = nvpa
        else if (ivpa < 1.0) then
            departure_outside_grid = .true. 
            ivpa = 1.0
        end if

    end subroutine calculate_zed_vpa_departure_idx



    subroutine advance_parallel_dynamics (pdf, phi)
        use finite_differences, only: second_order_centered
        use zgrid, only: nzgrid, ntubes
        use zgrid, only: delzed
        use vpamu_grids, only: nvpa, vpa, vpa_max, dvpa
        use vpamu_grids, only: maxwell_vpa, maxwell_mu
        use species, only: spec
        use stella_time, only: code_dt
        use stella_layouts, only: vmu_lo, kymus_lo
        use stella_layouts, only: iky_idx
        use redistribute, only: scatter, gather
        use dist_redistribute, only: kymus2vmus
        use arrays_dist_fn, only: g_kymus
        use extended_zgrid, only: nsegments, nzed_segment, neigen
        use extended_zgrid, only: map_to_extended_zgrid, map_from_extended_zgrid
        use extended_zgrid, only: map_to_iz_ikx_from_izext
        use fields, only: advance_fields
        ! use array_fields, only: apar, bpar (this is throwing an error)
        ! call advance_fields(gnew, phi, apar, bpar, dist = 'g')
    
        implicit none

        complex, dimension (:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (inout) :: pdf
        complex, dimension (:, :, -nzgrid:, :), intent (inout) :: phi

        integer :: ikymus, iky, imu, is
        integer :: ikx, iz, it, iv, ie
        integer :: izext, nz_ext
        integer :: ulim   ! dummy variable
        integer :: ia = 1 ! does not support flux annulus
        integer, dimension (:), allocatable :: iz_from_izext, ikx_from_izext
        complex, dimension (:, :), allocatable :: g_ext_1, g_ext_2
        complex, dimension (:), allocatable :: phi_ext, dphi_dz

        ! Used for interpolating g
        real :: dep_izext, dep_iv
        real, dimension(2) :: point
        integer :: grid_v_start, grid_zext_start, l, m 
        complex, dimension(4, 4) :: local4x4
        
    
        ! the input pdf is in the vmu_lo; re-map to work in kymus_lo
        call scatter(kymus2vmus, pdf, g_kymus)

        ! update g = <delta f> to account for motion along the lowest-order characteristic,
        ! for which the particle kinetic energy is constant
        do ikymus = kymus_lo%llim_proc, kymus_lo%ulim_proc
            ! map to the extended zed domain to ease calculations
            iky = iky_idx(kymus_lo, ikymus)
            do it = 1, ntubes
                do ie = 1, neigen(iky)
                    ! nz_ext is the number of grid points in the extended zed domain
                    nz_ext = nsegments(ie, iky) * nzed_segment + 1
                    allocate (iz_from_izext(nz_ext))
                    allocate (ikx_from_izext(nz_ext))
                    ! g_ext_1 and g_ext_2 will contain slices of the pdf on the extended zed domain
                    allocate (g_ext_1(nz_ext, nvpa)) 
                    allocate (g_ext_2(nz_ext, nvpa))

                    do iv = 1, nvpa
                        ! map from (kx, z, tube) to the extended zed domain; NB: ulim is a dummy argument
                        call map_to_extended_zgrid(it, ie, iky, g_kymus(:, :, :, iv, ikymus), g_ext_1(:, iv), ulim)
                    end do

                    call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)
                    ! set the pdf(t+dt) at each (zext, vpa) grid point to be equal to
                    ! the value of the pdf(t) at the (zext, vpa) that connects to it along the
                    ! characteristic with constant particle kinetic energy
                    ! NB: currently using a crude nearest-neighbour interpolation for departure point
                    do iv = 1, nvpa
                        do izext = 1, nz_ext
                            iz = iz_from_izext(izext)
                            ikx = ikx_from_izext(izext)

                            if (departure_point_outside_grid(ikx, iz, it, iv, ikymus)) then
                                g_ext_2(izext, iv) = 0.
                            else
                                ! Obtain interpolated values of g from g_ext_1 using departure izext and iv
                                dep_izext = departure_point_izext(ikx, iz, it, iv, ikymus)
                                dep_iv = departure_point_iv(ikx, iz, it, iv, ikymus)

                                ! Handle edge cases: top, bottom, left, and right points (see below)

                                ! Normal case:      Edge case:
                                ! o o o o o o o o   o o 1 x x 1 o o o 
                                ! o 1 1 1 1 o o o   o o 1 x x 1 o o o
                                ! o 1 x x 1 o o o   o o 1 1 1 1 o o o  
                                ! o 1 x x 1 o o o   o o 1 1 1 1 o o o 
                                ! o 1 1 1 1 o o o   o o o o o o o o o
                                ! o o o o o o o o   o o o o o o o o o 
                                
                                !> CHECK LOGIC AGAIN TO BE SURE !!

                                if (int(dep_iv) == 1) then
                                    grid_v_start = 1
                                    point(2) = modulo(dep_iv, 1.0) - 1.0
                                else if (int(dep_iv) == nvpa - 1) then
                                    grid_v_start = nvpa - 3
                                    point(2) = modulo(dep_iv, 1.0) + 1.0
                                else
                                    grid_v_start = int(dep_iv) - 1
                                    point(2) = modulo(dep_iv, 1.0)
                                end if

                                if (int(dep_izext) == 1) then
                                    grid_zext_start = 1
                                    point(1) = modulo(dep_izext, 1.0) - 1.0
                                else if (int(dep_izext) == nz_ext - 1) then
                                    grid_zext_start = nz_ext - 3
                                    point(1) = modulo(dep_izext, 1.0) + 1.0
                                else
                                    grid_zext_start = int(dep_izext) - 1
                                    point(1) = modulo(dep_izext, 1.0)
                                end if
                                
                                ! Construct grid array for interpolation
                                do m = 1, 4
                                    do l = 1, 4
                                        local4x4(l,m) = g_ext_1(grid_zext_start + l - 1, grid_v_start + m - 1)
                                    end do 
                                end do
                            
                                g_ext_2(izext, iv) = bicubic_spline(local4x4, point)

                            end if
                        end do
                    end do
                    
                    ! (Can I assume electrostatic field?) If so, the acceleration due to background Maxwellian population does not
                    ! change the population density of species at each spatial location. Therefore, any changes to phi is entirely
                    ! due to advection which we have already calculated previously. Therefore, use the midpoint to advance g due to    
                    ! this acceleration term? (i.e. we can compute dphi/dz in both the previous time-step as well as the current time-step)
                    ! Technically, we can even split this calculation into sub-timesteps
                    ! (i.e. advect -> dphi/dt -> advect -> dphi/dt -> etc. or maybe swap orders if that's more accurate), but 
                    ! this might defeat the whole purpose of using semi-lagrange (although there's still more freedom to choose larger substeps)
                    
                    ! call advance_fields(gnew, phi, apar, bpar, dist = 'g')
                    
                    ! phi_ext will contain a slice of the electrostatic potential on the extended zed domain
                    allocate (phi_ext(nz_ext))
                    ! dphi_dz will contain a slice of dphi/dz on the extended zed domain
                    allocate (dphi_dz(nz_ext))
                    ! map from (kx, z, tube) to the extended zed domain; NB: ulim is a dummy argument
                    call map_to_extended_zgrid(it, ie, iky, phi(iky, :, :, :), phi_ext, ulim)
                    ! compute dphi/dz using centered differences, with zero BCs
                    call second_order_centered(1, phi_ext, delzed(0), dphi_dz)
                    ! calculate the change in g over time dt due to the parallel acceleration of
                    ! particles in the background Maxwellian population (by the parallel electric field);
                    ! for now, use the parallel electric field at the departure point to obtain
                    ! the acceleration
                    do iv = 1, nvpa
                        do izext = 1, nz_ext
                            iz = iz_from_izext(izext)
                            ikx = ikx_from_izext(izext)
                            !> FLAG: This is WRONG, but just wanted to compile this, so rounded izext and iv
                            g_ext_2(izext, iv) = g_ext_2(izext, iv) + code_dt * spec(is)%zstm &
                            * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) &
                            * (- vpa_max + dvpa * (departure_point_iv(ikx, iz, it, iv, ikymus) - 1)) &
                            * dphi_dz(int(departure_point_izext(ikx, iz, it, iv, ikymus)))
                        end do
                    end do
                    do iv = 1, nvpa
                        ! map from the extended zed domain to (kx, z, tube); NB: ulim is a dummy argument
                        call map_from_extended_zgrid(it, ie, iky, g_ext_2(:, iv), g_kymus(:, :, :, iv, ikymus))
                    end do
                    ! deallocate g_ext so that it can be re-allocated with different size
                    deallocate (g_ext_1, g_ext_2, phi_ext, dphi_dz)
                    deallocate (iz_from_izext, ikx_from_izext)
                end do
            end do
        end do
    
        ! the output pdf should be in the vmu_lo; re-map from kymus_lo to vmu_lo
        call gather(kymus2vmus, g_kymus, pdf)

    end subroutine advance_parallel_dynamics

    subroutine finish_parallel_dynamics
        implicit none
        if (allocated(departure_point_izext)) deallocate(departure_point_izext)
        if (allocated(departure_point_iv)) deallocate(departure_point_iv)
        if (allocated(departure_point_outside_grid)) deallocate(departure_point_outside_grid)
    end subroutine finish_parallel_dynamics




    ! ##### Interpolating utilities #####!
    function bicubic_spline(grid_points, point) result(value)
        ! This function takes in a 4x4 array (grid_points) and point : [x,y] such that x and y are between -1 and 2, and returns
        ! a bicubic spline interpolation value for [x, y]

        ! Assumed uniform spacing, grid point coordinates are normalised to -1, 0, 1, 2 
        ! Not-a-knot spline used (seems to be default for spline interpolation in most packages)
        ! For simplicity, only a 4x4 grid has been splined here, which is valid if the function is sufficiently smooth over the whole grid
        ! This results in the nice property that we only have one fitted cubic polynomial in each axis, so we can reuse the same 
        ! expression across the whole region

        implicit none
        
        complex, intent(in), dimension(4, 4) :: grid_points
        real, intent(in), dimension(2) :: point
        real, dimension(4,4) :: m
        real :: real_value, imag_value
        complex :: value
        integer :: i, j, k, l

        real_value = 0.0
        imag_value = 0.0

        m = reshape([0.0, -1.0/3.0,  1.0/2.0,  -1.0/6.0, &
                     1.0, -1.0/2.0,   -1.0,     1.0/2.0, &
                     0.0,   1.0,     1.0/2.0,  -1.0/2.0, &
                     0.0, -1.0/6.0,    0.0,     1.0/6.0  ], [4,4])

        do i = 1, 4
            do j = 1, 4
                do k = 1, 4
                    do l = 1, 4
                        real_value = real_value + point(2)**(i-1) * m(i, j) * point(1)**(l-1) * m(l, k) * real(grid_points(k, j))
                        imag_value = imag_value + point(2)**(i-1) * m(i, j) * point(1)**(l-1) * m(l, k) * aimag(grid_points(k, j))
                    end do
                end do
            end do
        end do

        value = cmplx(real_value, imag_value)
    
    end function bicubic_spline

    function linear_interpolate(points, x) result(value)
        implicit none
        real, intent(in), dimension(2) :: points
        real, intent(in) :: x

        real :: value

        value = points(1) + (points(2)-points(1))*x

    end function linear_interpolate





end module parallel_dynamics
