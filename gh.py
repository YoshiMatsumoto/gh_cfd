import Rhino
import scriptcontext
import System.Guid
import rhinoscriptsyntax as rs


class CFD:

    def __init__(self, WX, WY, delta_t, Re, omega):
        self.WX = WX
        self.WY = WY
        self.vx = [[0] * (self.WX+1) for i in range(self.WY)]
        self.vx_after = [[0] * (self.WX+1) for i in range(self.WY)]
        self.vy = [[0] * self.WX for i in range(self.WY+1)]
        self.vy_after = [[0] * self.WX for i in range(self.WY+1)]
        self.s = [[0] * self.WX for i in range(self.WY)]

        self.p = [[0] * self.WX for i in range(self.WY)]
        self.p_after = [[0] * self.WX for i in range(self.WY)]
        self.rys =  [[0] * 1024 for i in range(2)]
        self.delta_t = delta_t                                  #デルタT
        self.Re = Re                                   #レイノルズ数
        self.omega = omega                                      #加速係数

    #移流
    def Adve(self):
        #vx の更新
        for i in range(self.WX-2):
            for j in range(self.WY-2):
                u = self.vx[i][j]
                v = (self.vx[i - 1][j] + self.vx[i][j] + self.vx[i - 1][j + 1] + self.vx[i][j + 1]) / 4
                if u>= 0 and v >= 0:
                    self.vx_after[i][j] = self.vx[i][j] - u*(self.vx[i][j] - self.vx[i -1][j])*self.delta_t - v * (self.vx[i][j] - self.vx[i][j-1]) * self.delta_t
                elif u < 0 and v >= 0:
                    self.vx_after[i][j] = self.vx[i][j] - u*(self.vx[i + 1][j] - self.vx[i][j])*self.delta_t - v * (self.vx[i][j] - self.vx[i][j-1]) * self.delta_t
                elif u >= 0 and v < 0:
                    self.vx_after[i][j] = self.vx[i][j] - u*(self.vx[i][j] - self.vx[i - 1][j])*self.delta_t - v * (self.vx[i][j + 1] - self.vx[i][j]) * self.delta_t
                elif u < 0 and v < 0:
                    self.vx_after[i][j] = self.vx[i][j] - u*(self.vx[i + 1][j] - self.vx[i][j])*self.delta_t - v * (self.vx[i][j + 1] - self.vx[i][j]) * self.delta_t

        #vy の更新
        for i in range(self.WX-2):
            for j in range(self.WY-2):
                u = (self.vy[i][j - 1] + self.vy[i + 1][j - 1] + self.vy[i][j] + self.vy[i + 1][j]) / 4
                v = self.vy[i][j]

                if u>= 0 and v >= 0:
                    self.vy_after[i][j] = self.vy[i][j] - u*(self.vy[i][j] - self.vy[i -1][j])*self.delta_t - v * (self.vy[i][j] - self.vy[i][j-1]) * self.delta_t
                elif u < 0 and v >= 0:
                    self.vy_after[i][j] = self.vy[i][j] - u*(self.vy[i + 1][j] - self.vy[i][j])*self.delta_t - v * (self.vy[i][j] - self.vy[i][j-1]) * self.delta_t
                elif u >= 0 and v < 0:
                    self.vy_after[i][j] = self.vy[i][j] - u*(self.vy[i][j] - self.vy[i - 1][j])*self.delta_t - v * (self.vy[i][j + 1] - self.vy[i][j]) * self.delta_t
                elif u < 0 and v < 0:
                    self.vy_after[i][j] = self.vy[i][j] - u*(self.vy[i + 1][j] - self.vy[i][j])*self.delta_t - v * (self.vy[i][j + 1] - self.vy[i][j]) * self.delta_t

        return

    #粘性
    def Viscosity(self):
        for i in range(self.WX-2):
            for  j in range(self.WY-2):
                self.vx_after[i][j] = self.vx[i][j] -1 / self.Re * (self.vx[i + 1][j] + self.vx[i][j + 1] + self.vx[i - 1][j] + self.vx[i][j - 1])*self.delta_t
                self.vy_after[i][j] = self.vy[i][j] -1 / self.Re * (self.vy[i + 1][j] + self.vy[i][j + 1] + self.vy[i - 1][j] + self.vy[i][j - 1])*self.delta_t
        return

    # 壁の設定
    def Set(self):
        for i in range(self.WX - 1):
            for j in range(self.WY - 1):
                if i == 0 or i == self.WX-1 or j == 0 or j == self.WY - 1:
                    self.vx[i][j] = 0
                    self.vx[i + 1][j] = 0
                    self.vy[i][j] = 0
                    self.vy[i][j + 1] = 0
        return

    # ダイバージェンスの計算
    def Div(self):
        for i in range(self.WX-2):
            for j in range(self.WY-2):
                self.s[i][j] = (-self.vx[i][j] - self.vy[i][j] + self.vx[i+1][j] + self.vy[i][j+1]) / self.delta_t
        return

    # 圧力のポアソン方程式
    def Poisson(self):
        for n in range(100):
            for i in range(self.WX-2):
                for j in range(self.WY-2):
                    if i == 1:
                        self.p[i - 1][j] = self.p[i][j]
                    elif i == self.WX-2:
                        self.p[i + 1][j] = self.p[i][j]
                    elif j == 1:
                        self.p[i][j - 1] = self.p[i][j]
                    elif j == self.WY-2:
                        self.p[i][j + 1] = self.p[i][j]

                    self.p[i][j] = (1 - self.omega)*self.p[i][j] + self.omega/4*(self.p[i - 1][j] + self.p[i + 1][j] + self.p[i][j - 1] + self.p[i][j + 1] - self.s[i][j])
        return

    #修正項
    def Rhs(self):
        for i in range(self.WX-2):
            for j in range(self.WY-2):
                self.vx[i][j] -= (self.p[i][j] - self.p[i - 1][j]) * self.delta_t
                self.vy[i][j] -= (self.p[i][j] - self.p[i][j - 1]) * self.delta_t
        return

# メッシュの作成
# def AddMesh():
#     mesh = Rhino.Geometry.Mesh()
#     # メッシュの頂点
#     mesh.Vertices.Add(0.0, 0.0, 1.0) #0
#     mesh.Vertices.Add(1.0, 0.0, 1.0) #1
#     mesh.Vertices.Add(2.0, 0.0, 1.0) #2
#     mesh.Vertices.Add(3.0, 0.0, 0.0) #3
#     mesh.Vertices.Add(0.0, 1.0, 1.0) #4
#     mesh.Vertices.Add(1.0, 1.0, 2.0) #5
#     mesh.Vertices.Add(2.0, 1.0, 1.0) #6
#     mesh.Vertices.Add(3.0, 1.0, 0.0) #7
#     mesh.Vertices.Add(0.0, 2.0, 1.0) #8
#     mesh.Vertices.Add(1.0, 2.0, 1.0) #9
#     mesh.Vertices.Add(2.0, 2.0, 1.0) #10
#     mesh.Vertices.Add(3.0, 2.0, 1.0) #11
#     # メッシュのフェイスの作成
#     mesh.Faces.AddFace(0, 1, 5, 4)
#     mesh.Faces.AddFace(1, 2, 6, 5)
#     mesh.Faces.AddFace(2, 3, 7, 6)
#     mesh.Faces.AddFace(4, 5, 9, 8)
#     mesh.Faces.AddFace(5, 6, 10, 9)
#     mesh.Faces.AddFace(6, 7, 11, 10)
#     # メッシュの法線を作成
#     mesh.Normals.ComputeNormals()
#     # メッシュの結合
#     mesh.Compact()

#     return mesh

# サーフェスを取得し、その範囲をシミュ範囲とする
# rec：サーフェス範囲
# div：分割単位（単位はRhinocerosのスケールに準拠）
domainU = int(round(rs.SurfaceDomain(rec, 0)[1]/div))
domainV = int(round(rs.SurfaceDomain(rec, 1)[1]/div))

# コンストラクタをシミュ範囲で定義する
# dt：ステップ間隔
# Re：レイノルズ定数
# omega；加速度
cfd = CFD(domainU, domainV, dt, Re, omega)

# step回シミュレーションを回す
for i in range(1024):
    cfd.Adve()
    cfd.Viscosity()
    cfd.Set()
    cfd.Div()
    cfd.Poisson()
    cfd.Rhs()