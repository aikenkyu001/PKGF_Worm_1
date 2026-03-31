import numpy as np
import math

# --- システム定数 (Fortran pkgf_constants と完全同期) ---
DIM = 32
N_NODES = 302
DT = 0.02
PI = 3.14159265358979323846
ALPHA = 1.0 / float(DIM)
VISCOSITY = 5.0

def invert_matrix(A):
    """DIMxDIM 行列の逆行列計算 (Fortran pkgf_math_utils.invert_matrix と同一アルゴリズム)"""
    n = A.shape[0]
    tA = A.copy()
    Ainv = np.eye(n, dtype=np.float64)
    for i in range(n):
        p = i
        pivot = abs(tA[i, i])
        for j in range(i + 1, n):
            if abs(tA[j, i]) > pivot:
                p = j
                pivot = abs(tA[j, i])
        if p != i:
            tA[[i, p]] = tA[[p, i]]
            Ainv[[i, p]] = Ainv[[p, i]]
        f = tA[i, i]
        tA[i, :] /= f
        Ainv[i, :] /= f
        for j in range(n):
            if i != j:
                f = tA[j, i]
                tA[j, :] -= f * tA[i, :]
                Ainv[j, :] -= f * Ainv[i, :]
    return Ainv

def get_determinant(A):
    """行列式の計算 (Fortran pkgf_math_utils.determinant と同一アルゴリズム)"""
    n = A.shape[0]
    tA = A.copy()
    det = 1.0
    for i in range(n):
        if abs(tA[i, i]) < 1.0e-30:
            return 0.0
        det *= tA[i, i]
        for j in range(i + 1, n):
            f = tA[j, i] / tA[i, i]
            tA[j, i:] -= f * tA[i, i:]
    return det

def matrix_exp(A):
    """マトリックス指数関数 (Pade 近似: Fortran pkgf_math_utils.matrix_exp と同一)"""
    id_m = np.eye(DIM, dtype=np.float64)
    A2 = A @ A
    A4 = A2 @ A2
    A6 = A4 @ A2
    U = A @ (A6 + 15120.0 * A4 + 604800.0 * A2 + 17643225600.0 * id_m) * (1.0 / 30240.0)
    V = 30240.0 * A6 + 3326400.0 * A4 + 86486400.0 * A2 + 17643225600.0 * id_m
    num = V + U
    den = V - U
    return invert_matrix(den) @ num

def adjoint_transform(K, Omega, dt):
    """随伴変換 (Fortran pkgf_math_utils.adjoint_transform と同一)"""
    H = matrix_exp(Omega * dt)
    H_inv = invert_matrix(H)
    return H @ K @ H_inv

class Food:
    def __init__(self, pos, intensity):
        self.pos = np.array(pos, dtype=np.float64)
        self.intensity = float(intensity)
        self.consumed = False

    def consume(self, amount):
        if self.consumed:
            return 0.0
        actual = amount
        self.intensity -= actual
        if self.intensity <= 0.0:
            self.intensity = 0.0
            self.consumed = True
        return actual

class Trace:
    def __init__(self, pos):
        self.pos = np.array(pos, dtype=np.float64)
        self.age = 1.0

class EnvironmentalPKGF:
    def __init__(self):
        self.foods = [Food([0.0, 0.0], 0.0)] # 初期化
        self.foods[0].consumed = True
        self.traces = []
        self.n_traces = 0

    def add_food(self, pos, intensity):
        self.foods[0] = Food(pos, intensity)

    def add_trace(self, pos):
        if len(self.traces) < 200:
            self.traces.append(Trace(pos))
        else:
            self.traces.pop(0)
            self.traces.append(Trace(pos))
        self.n_traces = len(self.traces)

    def get_curvature(self, obs_p):
        F = np.zeros((2, 2), dtype=np.float64)
        B_field = 0.0
        if not self.foods[0].consumed:
            diff = self.foods[0].pos - obs_p
            d_sq = np.sum(diff**2) + ALPHA**2
            B_field = self.foods[0].intensity / d_sq
        F[0, 1] = B_field
        F[1, 0] = -B_field
        return F

    def get_gradient(self, obs_p):
        grad = np.zeros(2, dtype=np.float64)
        if not self.foods[0].consumed:
            diff = self.foods[0].pos - obs_p
            d_sq = np.sum(diff**2) + ALPHA**2
            grad += (2.0 * self.foods[0].intensity / (d_sq**2)) * diff
        for trace in self.traces:
            diff = obs_p - trace.pos
            d_sq = np.sum(diff**2) + ALPHA**2
            grad += (2.0 * trace.age / (d_sq**2)) * diff
        return grad

    def step(self):
        # 減衰廃止 (Fortran と同期)
        pass

class PKGFWormNode:
    def __init__(self, node_id, segment_id, personality, nt, st, s_ph):
        self.id = node_id + 1
        self.segment_id = segment_id
        self.personality = personality
        self.nt = nt
        self.st = st
        self.pos = np.zeros(2, dtype=np.float64)
        self.v = np.zeros(DIM, dtype=np.float64)
        self.ext_F = np.zeros((DIM, DIM), dtype=np.float64)
        self.itf_nabla_s = np.zeros(DIM, dtype=np.float64)
        self.x = np.zeros(DIM, dtype=np.float64)
        
        # 初期状態 (Fortran pkgf_worm_node_mod.node_init と同一)
        for d in range(DIM):
            self.x[d] = 0.001 * (self.id + d + 1) / (N_NODES + DIM)
        self.x[14] = 0.1 * math.cos(s_ph) # Fortran x(15)
        self.x[15] = 0.1 * math.sin(s_ph) # Fortran x(16)
        
        # 並行鍵 K
        self.K = np.eye(DIM, dtype=np.float64)
        c_c, s_c = math.cos(PI/4.0), math.sin(PI/4.0)
        # Motor Coupling (Indices: 0-based for Fortran 1-based)
        for d in range(2):
            self.K[d, d] = c_c
            self.K[16+d, d] = -s_c
            self.K[d, 16+d] = s_c
            self.K[16+d, 16+d] = c_c
        
        # Knowledge R (x(17:24) sector is indices 16:24)
        # Fortran: self%K(17:24, 17:24) = matmul(self%K(17:24, 17:24), R) where R is Identity
        # R is identity here so no change needed.

    def hunger(self):
        # Fortran: sum(abs(self%x(25:32))) -> Python: x[24:32]
        return np.sum(np.abs(self.x[24:32])) / (DIM/4.0)

    def tension(self):
        # Fortran: sum(abs(self%x(1:8))) -> Python: x[0:8]
        return np.sum(np.abs(self.x[0:8])) / (DIM/4.0)

    def get_warp_factors(self):
        phi = math.sqrt(np.sum(self.x**2)) + math.sqrt(np.sum(self.itf_nabla_s**2))
        s = math.exp(phi * ALPHA)
        nabla_phi = self.itf_nabla_s * ALPHA
        nabla_phi[24:32] += ALPHA # Fortran 25:32
        return s, nabla_phi

    def get_v_pkgf(self):
        s, nabla_phi = self.get_warp_factors()
        phase_x = math.atan2(self.x[15], self.x[14]) # Fortran x(16), x(15)
        
        internal_F = np.zeros((DIM, DIM), dtype=np.float64)
        # Drive Manifold: Fortran (31, 17) -> Python (30, 16)
        internal_F[30, 16] = self.x[30] * ALPHA
        internal_F[16, 30] = -internal_F[30, 16]
        
        # Weathervane Steering: Fortran (17, 18) -> Python (16, 17)
        internal_F[16, 17] = math.sin(phase_x) * self.itf_nabla_s[16] - math.cos(phase_x) * self.itf_nabla_s[17]
        internal_F[17, 16] = -internal_F[16, 17]
        
        # Metabolism: Fortran (1, 31) -> Python (0, 30)
        internal_F[0, 30] = self.tension() * ALPHA
        internal_F[30, 0] = -internal_F[0, 30]
        
        # 自発的振動曲率: Fortran (15, 16) -> Python (14, 15)
        internal_F[14, 15] = self.tension() + ALPHA
        internal_F[15, 14] = -internal_F[14, 15]
        
        total_F = self.ext_F + internal_F
        Kx = self.K @ self.x
        v_target = np.zeros(DIM, dtype=np.float64)
        
        for d in range(DIM):
            # Viscosity applied to physical dimensions 17, 18 (Indices 16, 17)
            s_eff = s + (VISCOSITY if d in [16, 17] else 0.0)
            term_F = -np.dot(total_F[d, :], Kx)
            term_nabla = np.dot(self.K[d, :], nabla_phi)
            v_target[d] = (term_F + term_nabla) / s_eff
            
        return v_target

    def update(self, v_raw, dt):
        self.v = v_raw
        self.x += self.v * dt
        
        s, nabla_phi = self.get_warp_factors()
        Omega = np.zeros((DIM, DIM), dtype=np.float64)
        for i in range(DIM):
            for j in range(DIM):
                Omega[i, j] = (nabla_phi[i] * self.v[j] - nabla_phi[j] * self.v[i]) / s
        
        # Proprioceptive Feedback: Fortran (17, 1) -> Python (16, 0)
        v_mag = math.sqrt(np.sum(self.v[16:18]**2)) * ALPHA
        Omega[16, 0] += v_mag
        Omega[0, 16] -= v_mag
        
        self.K = adjoint_transform(self.K, Omega, dt)
        self.pos += self.v[16:18] * dt

class PKGFWormEngine:
    def __init__(self):
        self.nodes = []
        self.env = EnvironmentalPKGF()
        self.W = np.zeros((N_NODES, N_NODES), dtype=np.float64)
        
        # コネクトームと重みの初期化 (Fortran pkgf_connectome_data と同一)
        names = [f"N{i:03d}" for i in range(1, N_NODES + 1)]
        names[0] = "ALML"; names[1] = "ALMR"
        names[2] = "AVAL"; names[3] = "AVAR"
        names[4] = "AVBL"; names[5] = "AVBR"
        names[6] = "DB1";  names[7] = "VB1"
        
        for i in range(N_NODES):
            nt, st = "Interneuron", "Neutral"
            if names[i][:3] == "ALM":
                nt = "Sensory"; st = "Neutral"
            elif names[i][:2] in ["VB", "DB"]:
                nt = "Motor"; st = "Excitatory"
            elif names[i][:2] == "AV":
                nt = "Interneuron"; st = "Command"
            
            pos_z = float(i) / float(N_NODES)
            s_ph = 2.0 * PI * pos_z
            self.nodes.append(PKGFWormNode(i, i // 23, "Neutral", nt, st, s_ph))
            
        # Weights (Fortran indices W(j, i) are 1-based, Python 0-based)
        # ALML (1) connections
        self.W[0, 4] = 1.5; self.W[0, 2] = 0.5; self.W[0, 1] = 1.0 
        # ALMR (2) connections
        self.W[1, 5] = 1.5; self.W[1, 3] = 0.5; self.W[1, 0] = 1.0 
        # AVBL (5) -> VB1 (8)
        self.W[4, 7] = 2.0; self.W[4, 5] = 1.0 
        # AVBR (6) -> DB1 (7)
        self.W[5, 6] = 2.0; self.W[5, 4] = 1.0 
        # AVAL/AVAR
        self.W[2, 3] = 1.0; self.W[3, 2] = 1.0 
        
        for i in range(8, N_NODES):
            j = (i + 1) % N_NODES
            if j < 8: j = 8
            self.W[i, j] = 1.0

    def step(self, step_idx):
        if step_idx == 5:
            self.env.add_food([1.0, 1.0], 1.0)
            
        for i, node in enumerate(self.nodes):
            itf_raw = np.zeros(DIM, dtype=np.float64)
            for j in range(N_NODES):
                if i != j:
                    itf_raw += self.W[j, i] * (self.nodes[j].K @ self.nodes[j].x)
            node.itf_nabla_s = itf_raw * ALPHA
            
            if node.nt == "Sensory":
                node.ext_F = np.zeros((DIM, DIM), dtype=np.float64)
                F_env = self.env.get_curvature(node.pos)
                node.ext_F[16:18, 16:18] = F_env # Indices 16:18 for Fortran 17:18
                
                g_val = self.env.get_gradient(node.pos)
                # Fortran (31, 17) -> Python (30, 16)
                node.ext_F[30, 16] += g_val[0] * ALPHA
                node.ext_F[16, 30] = -node.ext_F[30, 16]
                node.ext_F[30, 17] += g_val[1] * ALPHA
                node.ext_F[17, 30] = -node.ext_F[30, 17]
                
                dist_sq = np.sum((node.pos - self.env.foods[0].pos)**2)
                feed_weight = math.exp(-dist_sq / (ALPHA**2))
                if not self.env.foods[0].consumed:
                    node.ext_F[30, 16] += feed_weight * ALPHA
                    node.ext_F[16, 30] = -node.ext_F[30, 16]
                    self.env.foods[0].consume(feed_weight * DT)
            
            v_raw = node.get_v_pkgf()
            node.update(v_raw, DT)
            
            if step_idx % 5 == 0 and i == 0:
                self.env.add_trace(node.pos)
        
        self.env.step()

if __name__ == "__main__":
    engine = PKGFWormEngine()
    print(" === PKGF-Worm Sync Verification (Python) ===")
    
    with open("../Logs/worm_log_python.txt", "w") as f:
        f.write("Step,HeadX,HeadY,Hunger,Tension,V,DetK\n")
        
        for step in range(301):
            n1 = engine.nodes[0]
            px, py = n1.pos[0], n1.pos[1]
            
            avg_h = sum(n.hunger() for n in engine.nodes) / N_NODES
            avg_t = sum(n.tension() for n in engine.nodes) / N_NODES
            avg_v = sum(math.sqrt(np.sum(n.v**2)) for n in engine.nodes) / N_NODES
            det_k = get_determinant(n1.K)
            
            # Fortran の書式に合わせる: (I5,A,F18.10,A,F18.10,A,F18.10,A,F18.10,A,F18.10,A,F18.10)
            f.write(f"{step:5d},{px:18.10f},{py:18.10f},{avg_h:18.10f},{avg_t:18.10f},{avg_v:18.10f},{det_k:18.10f}\n")
            f.flush() # 1ステップごとにディスクへ保存を確定
            
            if step % 50 == 0:
                print(f"Step {step:5d} | Pos: ({px:10.2f}, {py:10.2f}) | Hunger: {avg_h:10.6f} | V: {avg_v:10.6f}")
            
            engine.step(step)
            
    print(" Log saved to Logs/worm_log_python.txt")
