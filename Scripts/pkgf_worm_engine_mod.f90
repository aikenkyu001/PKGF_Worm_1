module pkgf_worm_engine_mod
    use pkgf_constants
    use pkgf_worm_node_mod
    use pkgf_environment_mod
    use pkgf_connectome_data
    implicit none
    
    type PKGFWormEngine
        type(PKGFWormNode) :: nodes(N_NODES)
        type(EnvironmentalPKGF) :: env
        real(DP) :: W(N_NODES, N_NODES)
    contains
        procedure :: init => engine_init
        procedure :: step => engine_step
    end type PKGFWormEngine
contains
    subroutine engine_init(self)
        class(PKGFWormEngine), intent(inout) :: self
        integer :: i
        real(DP) :: s_ph, pos_z
        character(len=20) :: nt, st
        print *, " [Engine] Deploying full 302-node connectome (Final Absolute Purity)..."
        call init_neuron_names()
        call init_weight_matrix(self%W)
        do i = 1, N_NODES
            nt = "Interneuron"; st = "Neutral"
            if (NEURON_NAMES(i)(1:3) == "ALM") then
                nt = "Sensory"; st = "Neutral"
            else if (NEURON_NAMES(i)(1:2) == "VB" .or. NEURON_NAMES(i)(1:2) == "DB") then
                nt = "Motor"; st = "Excitatory"
            else if (NEURON_NAMES(i)(1:2) == "AV") then
                nt = "Interneuron"; st = "Command"
            end if
            pos_z = real(i - 1, DP) / real(N_NODES, DP)
            s_ph = 2.0d0 * PI * pos_z
            call self%nodes(i)%init(i, i/23, "Neutral", nt, st, s_ph)
        end do
        self%env%foods(1)%consumed = .true.
        self%env%n_traces = 0
    end subroutine engine_init

    subroutine engine_step(self, step_idx)
        class(PKGFWormEngine), intent(inout) :: self
        integer, intent(in) :: step_idx
        integer :: i, j
        real(DP) :: itf_raw(DIM), v_raw(DIM)
        if (step_idx == 5) call self%env%add_food((/1.0d0, 1.0d0/), 1.0d0)
        do i = 1, N_NODES
            itf_raw = 0.0d0
            do j = 1, N_NODES
                if (i /= j) then
                    ! 他者の流束をありのままに集計。BETA による希釈を廃止。
                    itf_raw = itf_raw + self%W(j, i) * matmul(self%nodes(j)%K, self%nodes(j)%x)
                end if
            end do
            ! 干渉はALPHA（次元正規化）のみでスケール
            self%nodes(i)%itf_nabla_s = itf_raw * ALPHA
            
            if (self%nodes(i)%nt == "Sensory") then
                self%nodes(i)%ext_F = 0.0d0
                block
                    real(DP) :: F_env(2, 2), dist_sq, feed_weight, g_val(2)
                    F_env = self%env%get_curvature(self%nodes(i)%pos)
                    self%nodes(i)%ext_F(17:18, 17:18) = F_env
                    g_val = self%env%get_gradient(self%nodes(i)%pos)
                    self%nodes(i)%ext_F(31, 17) = self%nodes(i)%ext_F(31, 17) + g_val(1) * ALPHA
                    self%nodes(i)%ext_F(17, 31) = -self%nodes(i)%ext_F(31, 17)
                    self%nodes(i)%ext_F(31, 18) = self%nodes(i)%ext_F(31, 18) + g_val(2) * ALPHA
                    self%nodes(i)%ext_F(18, 31) = -self%nodes(i)%ext_F(31, 18)
                    dist_sq = sum((self%nodes(i)%pos - self%env%foods(1)%pos)**2)
                    feed_weight = exp(-dist_sq / (ALPHA**2)) 
                    if (.not. self%env%foods(1)%consumed) then
                        self%nodes(i)%ext_F(31, 17) = self%nodes(i)%ext_F(31, 17) + feed_weight * ALPHA
                        self%nodes(i)%ext_F(17, 31) = -self%nodes(i)%ext_F(31, 17)
                        block
                            real(DP) :: dummy_eat
                            dummy_eat = self%env%foods(1)%consume(feed_weight * DT) 
                        end block
                    end if
                end block
            end if
            call self%nodes(i)%get_v_pkgf(v_raw)
            call self%nodes(i)%update(v_raw, DT)
            if (mod(step_idx, 5) == 0 .and. i == 1) then
                call self%env%add_trace(self%nodes(i)%pos)
            end if
        end do
        call self%env%step()
    end subroutine engine_step
end module pkgf_worm_engine_mod
