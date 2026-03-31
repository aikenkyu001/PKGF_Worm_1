module pkgf_worm_node_mod
    use pkgf_constants
    use pkgf_math_utils
    implicit none
    
    type PKGFWormNode
        integer :: id, segment_id
        character(len=20) :: personality, nt, st
        real(DP) :: x(DIM), v(DIM), K(DIM, DIM)
        real(DP) :: pos(2)
        real(DP) :: ext_F(DIM, DIM) 
        real(DP) :: itf_nabla_s(DIM)
        contains
        procedure :: init => node_init
        procedure :: hunger => node_hunger
        procedure :: tension => node_tension
        procedure :: get_warp_factors => node_get_warp_factors
        procedure :: get_v_pkgf => node_get_v_pkgf
        procedure :: update => node_update
    end type PKGFWormNode

contains
    subroutine node_init(self, id, s_id, personality, nt, st, s_ph)
        class(PKGFWormNode), intent(inout) :: self
        integer, intent(in) :: id, s_id
        character(len=*), intent(in) :: personality, nt, st
        real(DP), intent(in) :: s_ph
        real(DP) :: c_c, s_c, R(8,8)
        integer :: d, i
        self%id = id; self%segment_id = s_id; self%personality = personality
        self%nt = nt; self%st = st
        self%pos = 0.0d0; self%ext_F = 0.0d0; self%itf_nabla_s = 0.0d0
        
        ! 初期状態: 以前の s_ph (外部位相) を x(15:16) の回転初期値として「多様体内に」埋め込む
        do d = 1, DIM
            self%x(d) = 0.001d0 * real(id + d, DP) / real(N_NODES + DIM, DP)
        end do
        self%x(15) = 0.1d0 * cos(s_ph)
        self%x(16) = 0.1d0 * sin(s_ph)
        self%v = 0.0d0
        
        ! 並行鍵 K: 初期状態の知識（構造）
        self%K = 0.0d0; do d = 1, DIM; self%K(d, d) = 1.0d0; end do
        c_c = cos(PI/4.0d0); s_c = sin(PI/4.0d0)
        do d = 0, 1
            ! Motor Coupling
            self%K(1+d, 1+d) = c_c; self%K(17+d, 1+d) = -s_c
            self%K(1+d, 17+d) = s_c; self%K(17+d, 17+d) = c_c
        end do
        ! 知識 R を x(17:24) セクターへ統合
        R = 0.0d0; do i = 1, 8; R(i, i) = 1.0d0; end do
        self%K(17:24, 17:24) = matmul(self%K(17:24, 17:24), R)
    end subroutine node_init

    function node_hunger(self) result(h)
        class(PKGFWormNode), intent(in) :: self
        real(DP) :: h
        h = sum(abs(self%x(25:32))) / real(DIM/4, DP)
    end function node_hunger

    function node_tension(self) result(t)
        class(PKGFWormNode), intent(in) :: self
        real(DP) :: t
        t = sum(abs(self%x(1:8))) / real(DIM/4, DP)
    end function node_tension

    subroutine node_get_warp_factors(self, s, nabla_phi)
        class(PKGFWormNode), intent(in) :: self
        real(DP), intent(out) :: s, nabla_phi(DIM)
        real(DP) :: phi
        ! 【真実】ポテンシャル phi は内的状態 x のエネルギー密度そのもの。
        phi = sqrt(sum(self%x**2)) + sqrt(sum(self%itf_nabla_s**2))
        s = exp(phi * ALPHA) 
        nabla_phi = self%itf_nabla_s * ALPHA
        nabla_phi(25:32) = nabla_phi(25:32) + ALPHA
    end subroutine node_get_warp_factors

    subroutine node_get_v_pkgf(self, v_target)
        class(PKGFWormNode), intent(in) :: self
        real(DP), intent(out) :: v_target(DIM)
        real(DP) :: s, nabla_phi(DIM), internal_F(DIM, DIM), total_F(DIM, DIM), phase_x
        integer :: d
        call self%get_warp_factors(s, nabla_phi)
        
        ! 【真実】位相は多様体座標 x(15), x(16) から動的に抽出される。
        phase_x = atan2(self%x(16), self%x(15))
        
        internal_F = 0.0d0
        ! Drive Manifold
        internal_F(31, 17) = self%x(31) * ALPHA
        internal_F(17, 31) = -internal_F(31, 17)
        
        ! Weathervane Steering: 振動位相 (phase_x) と勾配成分の純粋な幾何学的結合
        internal_F(17, 18) = sin(phase_x) * self%itf_nabla_s(17) - cos(phase_x) * self%itf_nabla_s(18)
        internal_F(18, 17) = -internal_F(17, 18)
        
        ! Metabolism
        internal_F(1, 31) = self%tension() * ALPHA
        internal_F(31, 1) = -internal_F(1, 31)
        
        ! 自発的振動曲率: 15-16平面における回転流動を生成（リズムの自発的創発）
        internal_F(15, 16) = self%tension() + ALPHA
        internal_F(16, 15) = -internal_F(15, 16)
        
        total_F = self%ext_F + internal_F
        
        do d = 1, DIM
            ! 全ノード、全次元で一律の等方的計量。動きの差は K の構造のみから生まれる。
            v_target(d) = (- dot_product(total_F(d, :), matmul(self%K, self%x)) + &
                             dot_product(self%K(d, :), nabla_phi)) / (s + merge(VISCOSITY, 0.0d0, d==17 .or. d==18))
        end do
    end subroutine node_get_v_pkgf

    subroutine node_update(self, v_raw, dt)
        class(PKGFWormNode), intent(inout) :: self
        real(DP), intent(in) :: v_raw(DIM), dt
        self%v = v_raw
        
        ! 【絶対真理】唯一の更新式: dx/dt = v
        self%x = self%x + self%v * dt
        
        block
            real(DP) :: Omega(DIM, DIM), ss, ns(DIM)
            integer :: i, j
            call self%get_warp_factors(ss, ns)
            Omega = 0.0d0
            do i = 1, DIM
                do j = 1, DIM
                    Omega(i, j) = (ns(i) * self%v(j) - ns(j) * self%v(i)) / ss
                end do
            end do
            ! Proprioceptive Feedback: 移動速度の内的緊張への還流
            Omega(17, 1) = Omega(17, 1) + sqrt(sum(self%v(17:18)**2)) * ALPHA
            Omega(1, 17) = -Omega(17, 1)
            
            call adjoint_transform(self%K, Omega, dt)
        end block
        self%pos = self%pos + self%v(17:18) * dt
        ! 独立変数 sp_p の更新を廃止。全ては x に集約された。
    end subroutine node_update
end module pkgf_worm_node_mod
