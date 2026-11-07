# mod object can be defined

    Code
      print(m)
    Output
      <mod object>
      
      Formula: choice ~ price | 1 | time
      Random effects: ASC
      Alternative names: A, B, C
      Coefficients: 4 fixed, 2 random
      Free parameters total: 10
      
      Utility equations:
      
        U_A = b[1] * price_A + b[2] * time_A + eps[1]
        U_B = b[1] * price_B + b[3] * time_B + b[5] + eps[2]
        U_C = b[1] * price_C + b[4] * time_C + b[6] + eps[3]
      
      Parameter restrictions:
      
        b = H_b theta_b
        length(theta_b): 6
        H_b (6 x 6) =
                  theta_1 theta_2 theta_3 theta_4 theta_5 theta_6
      b[1]:price        1       0       0       0       0       0
      b[2]:time_A       0       1       0       0       0       0
      b[3]:time_B       0       0       1       0       0       0
      b[4]:time_C       0       0       0       1       0       0
      b[5]:ASC_B        0       0       0       0       1       0
      b[6]:ASC_C        0       0       0       0       0       1
      
      
        O = vech(chol(Omega)) = H_O theta_O
        length(theta_O): 2
        H_O (3 x 2) =
                    theta_7 theta_8
      ASC_B_:_ASC_B       1       0
      ASC_C_:_ASC_B       0       0
      ASC_C_:_ASC_C       0       1
      
      
        L = vech(chol(Sigma)) = H_L theta_L + f_L
        length(theta_L): 2
        H_L (6 x 2) =
            theta_9 theta_10
      A_:_A       0        0
      B_:_A       0        0
      C_:_A       0        0
      B_:_B       0        0
      C_:_B       1        0
      C_:_C       0        1
        f_L (6 x 1) =
             
      A_:_A 0
      B_:_A 0
      C_:_A 0
      B_:_B 1
      C_:_B 0
      C_:_C 0
      

---

    Code
      m$print()
    Output
      <mod object>
      
      Formula: choice ~ price | 1 | time
      Random effects: ASC
      Alternative names: A, B, C
      Coefficients: 4 fixed, 2 random
      Free parameters total: 10
      
      Utility equations:
      
        U_A = b[1] * price_A + b[2] * time_A + eps[1]
        U_B = b[1] * price_B + b[3] * time_B + b[5] + eps[2]
        U_C = b[1] * price_C + b[4] * time_C + b[6] + eps[3]
      
      Parameter restrictions:
      
        b = H_b theta_b
        length(theta_b): 6
        H_b (6 x 6) =
                  theta_1 theta_2 theta_3 theta_4 theta_5 theta_6
      b[1]:price        1       0       0       0       0       0
      b[2]:time_A       0       1       0       0       0       0
      b[3]:time_B       0       0       1       0       0       0
      b[4]:time_C       0       0       0       1       0       0
      b[5]:ASC_B        0       0       0       0       1       0
      b[6]:ASC_C        0       0       0       0       0       1
      
      
        O = vech(chol(Omega)) = H_O theta_O
        length(theta_O): 2
        H_O (3 x 2) =
                    theta_7 theta_8
      ASC_B_:_ASC_B       1       0
      ASC_C_:_ASC_B       0       0
      ASC_C_:_ASC_C       0       1
      
      
        L = vech(chol(Sigma)) = H_L theta_L + f_L
        length(theta_L): 2
        H_L (6 x 2) =
            theta_9 theta_10
      A_:_A       0        0
      B_:_A       0        0
      C_:_A       0        0
      B_:_B       0        0
      C_:_B       1        0
      C_:_C       0        1
        f_L (6 x 1) =
             
      A_:_A 0
      B_:_A 0
      C_:_A 0
      B_:_B 1
      C_:_B 0
      C_:_C 0
      

