#ввод-вывод с использованием BioPython из текстового файла
#ООП
#тесты
#выводим всю таблицу динамического программирования
#последовательности длиной >10

#var: 1.4 Парная СММ, алгоритм Витерби, логарифмические счета
# [Vm, Vx, Vy]


from Bio import SeqIO
#from Bio import Alphabet
import numpy as np

class Alignment:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.opt = [] #сюда будут сохраняться выравненные последовательности
        self.seq = '' #последовательность, которая будет выравниваться на вторую
        self.len1 = len(seq1) + 1  # иначе не посчитает последнюю букву
        self.len2 = len(seq2) + 1
        if (len(seq1) < len(seq2)):
            for i in range(self.len2 - 1):
                self.opt.append([seq2[i], 0])
            self.seq = self.seq1
        else:
            for i in range(self.len1 - 1):
                self.opt.append([seq1[i], 0])
            self.seq = self.seq2

        #вероятности
        self.probability = {'AA': np.log(0.5), 'CC': np.log(0.5), 'GG': np.log(0.5), 'TT': np.log(0.5),
                       'CT': np.log(0.05), 'TC': np.log(0.05), 'AG': np.log(0.05), 'GA': np.log(0.05),
                       'AT': np.log(0.3), 'TA': np.log(0.3), 'GC': np.log(0.3), 'CG': np.log(0.3),
                       'GT': np.log(0.15), 'TG': np.log(0.15), 'AC': np.log(0.15), 'CA': np.log(0.15),
                       'q': np.log(0.25),
                       'delta': np.log(0.2),
                       'tau': np.log(0.1),
                       'eps': np.log(0.1),
                       'XM_YM': np.log(0.8), 'BM_MM': np.log(0.5)}

    def opt_alignment(self, matrix_dp):
        # ищем оптимальное выравнивание по max
        i = self.len1 - 1
        j = self.len2 - 1
        max_p = max(matrix_dp[i][j])
        k = len(self.seq) - 1  # индекс для обратного отсчета в меньшей послед.
        for l in range(len(self.opt) - 1, -1, -1):
            if max_p == matrix_dp[i][j][1]:  # Vx = max
                j -= 1
                self.opt[l][1] = '-'
            elif max_p == matrix_dp[i][j][2]:  # Vy = max
                self.opt[l][1] = '-'
                i -= 1
            else:
                self.opt[l][1] = self.seq[k]  # Vm = max или одинаковые
                k -= 1
                i -= 1
                j -= 1
            max_p = max(matrix_dp[i][j])

        print("Оптимальное выравнивание:")
        for i in range(len(self.opt)):
            print(self.opt[i][0], end=' ')
        print("\n")
        for i in range(len(self.opt)):
            print(self.opt[i][1], end=' ')

    def matrix_dp(self):

        matrix_dp = []
        for i in range(self.len1):
            matrix_dp.append([])
            for j in range(self.len2):
                matrix_dp[i].append([])

        # инициализация
        matrix_dp[0][0] = [0, -np.inf, -np.inf]

        matrix_dp[0][1].append(-np.inf)
        matrix_dp[0][1].append(round(matrix_dp[0][0][0] + self.probability['delta'] + self.probability['q'], 3))
        matrix_dp[0][1].append(-np.inf)

        matrix_dp[1][0].append(-np.inf)
        matrix_dp[1][0].append(-np.inf)
        matrix_dp[1][0].append(round(matrix_dp[0][0][0] + self.probability['delta'] + self.probability['q'], 3))

        for j in range(2, self.len2):
            matrix_dp[0][j].append(-np.inf)
            matrix_dp[0][j].append(round(matrix_dp[0][j - 1][1] + self.probability['eps'] + self.probability['q'], 3))
            matrix_dp[0][j].append(-np.inf)

        for i in range(2, self.len1):
            matrix_dp[i][0].append(-np.inf)
            matrix_dp[i][0].append(-np.inf)
            matrix_dp[i][0].append(round(matrix_dp[i - 1][0][2] + self.probability['eps'] + self.probability['q'], 3))


        for i in range(1, self.len1):
            for j in range(1, self.len2):
                k = 0

                # Vm
                Vm_max = matrix_dp[i - 1][j - 1][0] + self.probability['BM_MM']
                Vx_max = matrix_dp[i - 1][j - 1][1] + self.probability['XM_YM']
                Vy_max = matrix_dp[i - 1][j - 1][2] + self.probability['XM_YM']

                if (Vm_max >= Vx_max) and (Vm_max >= Vy_max):
                    k += Vm_max
                elif (Vx_max >= Vm_max) and (Vx_max >= Vy_max):
                    k += Vx_max
                else:
                    k += Vy_max

                k += self.probability[self.seq1[i - 1] + self.seq2[j - 1]]

                matrix_dp[i][j].append(round(k, 3))

                # Vx
                k = 0
                Vm_max = matrix_dp[i][j - 1][0] + self.probability['delta']
                Vx_max = matrix_dp[i][j - 1][1] + self.probability['eps']

                if (Vm_max >= Vx_max):
                    k += Vm_max + self.probability['q']
                else:
                    k += Vx_max + self.probability['q']

                matrix_dp[i][j].append(round(k, 3))

                # Vy
                k = 0
                Vm_max = matrix_dp[i - 1][j][0] + self.probability['delta']
                Vy_max = matrix_dp[i - 1][j][2] + self.probability['eps']

                if (Vm_max >= Vy_max):
                    k += Vm_max + self.probability['q']
                else:
                    k += Vy_max + self.probability['q']

                matrix_dp[i][j].append(round(k, 3))

        print("Матрица динамического программирования:\n")
        for l in range(self.len1):
            for k in range(self.len2):
                print(matrix_dp[l][k], end='  |  ')
            print('\n')

        print("Вес оптимального выравнивания:", round(max(matrix_dp[self.len1 - 1][self.len2 - 1]) + self.probability['tau'], 3)) #перепроверить, какая вероятность

        self.opt_alignment(matrix_dp)


def validate(seq):
    alphabet = 'ACGT'
    for i in seq:
        if i not in alphabet:
            return False
        else:
            return True


if __name__ == '__main__':
    file1 = open("2.fasta", "r")
    file2 = open("1.fasta", "r")

    seq1 = SeqIO.read(file1, "fasta")
    seq2 = SeqIO.read(file2, "fasta")

    val = validate(seq1.seq)
    if val == False:
        print("Последовательность содержит недопустимые буквы")
        exit()
    else:
        val = validate(seq2.seq)
        if val == False:
            print("Последовательность содержит недопустимые буквы")
            exit()
        else:
            align = Alignment(seq1.seq, seq2.seq)
            align.matrix_dp()



