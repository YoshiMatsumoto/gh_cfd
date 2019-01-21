class CFD:

    def __init__(self, WX, WY, delta_t, Re, omega):
        self.WX = WX
        self.WY = WY
        self.vx = [[0] * (self.WX+1) for i in range(self.WY)]
        self.vx_after = [[0] * (self.WX+1) for i in range(self.WY)]
        self.vy = [[0] * self.WX for i in range(self.WY+1)]
        self.vy_after = [[0] * self.WX for i in range(self.WY+1)]
        self.p = [[0] * self.WX for i in range(self.WY)]
        self.p_after = [[0] * self.WX for i in range(self.WY)]
        self.s = [[0] * self.WX for i in range(self.WY)]

        self.rys =  [[0] * 1024 for i in range(2)]
        self.delta_t = delta_t                                  #デルタT
        self.delta=x = delta_x
        self.delta_y = delta_y
        self.Re = Re                                   #レイノルズ数
        self.omega = omega                                      #加速係数
    
    #時間ステップ数の計算
    def deltaT(self, delta_x, delta_y, u, v):
        minx = 0.2 * (delta_x / u)
        miny = 0.2 * (delta_y / v)
        if minx >= miny:
            return minx
        else:
            return miny

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

def main():
    # コンストラクタ
    cfd = CFD(12, 18, 0.2, 10000, 1.8)

    # 1024回回す
    # for i in range(1024):

    #     # 風の発生地点
    #     cfd.vx[4][3] = 1
    #     cfd.Adve()
    #     cfd.Viscosity()
    #     cfd.Set()
    #     cfd.Div()
    #     cfd.Poisson()
    #     cfd.Rhs()
    # print(cfd.vx)
    # print(cfd.vy)


if __name__ == '__main__':
    main()