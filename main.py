# import Rhino
# import scriptcontext
# import System.Guid
# import rhinoscriptsyntax as rs

class CFD:
    def __init__(self):
        self.WX = 12                                        #壁x位置
        self.WY = 12                                        #壁y位置
        self.vx = [[0] * self.WX+1 for i in range(self.WY)]           #速度x成分[配列]
        self.vx_after = [[0] * self.WX+1 for i in range(self.WY)]     #修正後速度x成分
        self.vy = [[0] * self.WX for i in range(self.WY+1)]           #速度y成分
        self.vy_after = [[0] * self.WX for i in range(self.WY+1)]     #修正後速度y成分
        self.s = [[0] * self.WX for i in range(self.WY)]

        self.p = [[0] * self.WX for i in range(self.WY)]              #圧力
        self.p_after = [[0] * self.WX for i in range(self.WY)]        #修正後圧力
        self.rys =  [[0] * 1024 for i in range(2)]          #粒子のx,y座標
        self.delta_t = 0.2                                  #デルタT
        self.Re=1000000.0                                   #レイノルズ数
        self.omega=1.8                                      #加速係数

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
        for i in range(self.WX):
            for j in range(self.WY):
                if i == 0 or i == self.WX-1 or j == 0 or j == self.WY - 1:
                    self.vx[i][j] = 0
                    self.vx[i + 1][j] == 0
                    self.vy[i][j] == 0
                    self.vy[i][j + 1] == 0
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


def main():
    step = 1024

    while step == 0:
        step -= 1
        cfd = CFD()
        cfd.Adve()
        cfd.Viscosity()
        cfd.Set()
        cfd.Div()
        cfd.Poisson()
        cfd.Rhs()


if __name__ == '__main__':
    main()