from sortedcontainers import SortedDict
from matplotlib import pyplot as plt
import numpy as np


def statistical_bound_of_waiting_time(p, k, alpha=0.05):
    """

    :param p:  average accuracy of read
    :param k:
    :param alpha: statistical bound
    :return:statistical bound of waiting time
    """
    p_k = float(p) ** (k)
    qp_k = (1 - p) * p_k
    x = 0
    sum_0_xk1 = 0.00
    sum = 0

    last_k_prob = np.zeros((k + 1), dtype=np.double)

    while sum < 1 - alpha:
        if x < k:
            last_k_prob[x % (k + 1)] = 0.00
        elif x == k:
            last_k_prob[x % (k + 1)] = p_k
        else:
            sum_0_xk1 += last_k_prob[x % (k + 1)]
            last_k_prob[x % (k + 1)] = qp_k * (1 - sum_0_xk1)

        sum += last_k_prob[x % (k + 1)]
        x += 1

    return x


def find_floor_key(d, key):
    try:
        i_loc = d.index(key)
        if i_loc == 0:
            return None, None, None
        else:
            return i_loc - 1, d.iloc[i_loc - 1], d[d.iloc[i_loc - 1]]

    except ValueError:
        # just insert something
        d[key] = 0
        i_loc = d.index(key)
        if i_loc == 0:
            del d[key]
            return None, None, None
        else:
            del d[key]
            return i_loc - 1, d.iloc[i_loc - 1], d[d.iloc[i_loc - 1]]


def find_not_smaller_key(d, key):
    try:
        i_loc = d.index(key)
        return i_loc, d.iloc[i_loc], d[d.iloc[i_loc]]
    except ValueError:
        # just insert something
        d[key] = 0
        i_loc = d.index(key)
        if i_loc == len(d) - 1:
            del d[key]
            return None, None, None
        else:

            del d[key]
            return i_loc, d.iloc[i_loc], d[d.iloc[i_loc]]


def find_not_bigger_key(d, key):
    try:
        i_loc = d.index(key)
        return i_loc, d.iloc[i_loc], d[d.iloc[i_loc]]
    except ValueError:
        # just insert something
        d[key] = 0
        i_loc = d.index(key)
        if i_loc == 0:
            del d[key]
            return None, None, None
        else:

            del d[key]
            return i_loc - 1, d.iloc[i_loc - 1], d[d.iloc[i_loc - 1]]


class GroupHit(object):
    def __init__(self, group_line, group_size = 12):
        sp = group_line.strip().split("\t")
        self.target = sp[0]
        self.target_len = int(sp[1])
        # print(self.target_len)
        self.query = sp[2]
        self.query_len = int(sp[3])
        # print(self.query_len)
        self.aligned = False
        self.chain_align = None

        # parameter for output region
        self.query_ali_start = None
        self.query_ali_end = None
        self.target_ali_start = None
        self.target_ali_end = None

        self.aligned_base = 0
        self.overlap_length = 0
        self.query_overlap_start = None
        self.query_overlap_end = None
        self.target_overlap_start = None
        self.target_overlap_end = None
        self.group_size = group_size    # default group size to be considered as a valid group

        if int(sp[4]) == 0:
            self.forward = True

        else:
            self.forward = False

        self.ratio = 0.5

        # L and I is for DP chaining
        self.I = []

        groups = []
        for group in sp[5:]:
            hits = []
            group_sp = group.strip(",").split(",")
            if len(group_sp) > 1 or int(group_sp[0].split()[2]) >= self.group_size:
                for hit in group_sp:
                    hit_sp = hit.split(" ")

                    if self.forward:
                        hits.append((int(hit_sp[0]), self.target_len - int(hit_sp[1]) + int(hit_sp[0]), int(hit_sp[2])))

                    else:
                        hits.append((int(hit_sp[0]), int(hit_sp[1]) - int(hit_sp[0]), int(hit_sp[2])))

                self.I.append((hits[0][0] - hits[0][2], len(groups), 0))
                self.I.append((hits[-1][0], len(groups), -1))
                groups.append(hits)

        self.groups = groups

    def fw_chain_groups(self):

        r = len(self.I)
        # L = FastRBTree()
        L_dict = SortedDict()
        V = [0] * len(self.groups)
        back_track = [-1] * len(self.groups)

        self.I.sort()
        # print self.I
        ## print self.groups

        for i in range(r):
            # I_list[i] is a start point, noticed we go through the chain from botton to top
            if self.I[i][2] == 0:
                k = self.I[i][1]

                l_k = self.groups[k][0][1] - self.groups[k][0][2]
                start_y = l_k
                start_x = self.I[i][0]
                # print(start_x, start_y)
                end_x = self.groups[k][-1][0]
                end_y = self.groups[k][-1][1]
                diagonal = start_y - start_x
                # find largest h_j strictly smaller than l_k and also not off diagonal

                j_index, j_key, j_value = find_floor_key(L_dict, l_k)
                # print("j", j_value)
                while j_index != None:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = j_value[1]
                    j_x = self.groups[j][-1][0]
                    j_y = self.groups[j][-1][1]
                    d_x = start_x - j_x
                    d_y = start_y - j_y
                    # print("start", j_x, j_y)
                    # print(d_x, d_y)
                    if min(d_x, d_y) / max(d_x, d_y) > self.ratio:

                        prev_score = V[j]
                        break

                    else:
                        j_index, j_key, j_value = find_floor_key(L_dict, j_y)


                else:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = -1
                    prev_score = 0

                V[k] = prev_score + v_j
                back_track[k] = j

            # is a end point
            else:

                k = self.I[i][1]
                h_k = self.groups[k][-1][1]

                l_k = self.groups[k][0][1] - self.groups[k][0][2]
                start_y = l_k
                start_x = self.I[i][0]
                diagonal = start_y - start_x

                # print('start', start_x, start_y)

                j_index, j_key, j_value = find_not_smaller_key(L_dict, h_k)
                if j_index != None:
                    # print("1st")
                    j = j_value[1]

                    j_x = self.groups[j][-1][0]
                    j_y = self.groups[j][-1][1]
                    j_diagonal = j_y - j_x

                    V_j = j_value[0]

                    if V[k] > V_j or abs(j_diagonal - diagonal) > 10:
                        L_dict[h_k] = (V[k], k)

                else:
                    if len(L_dict) == 0:
                        L_dict[h_k] = (V[k], k)

                    else:
                        # print(L_dict.peekitem())
                        last_j = L_dict.peekitem()[1][1]
                        j_x = self.groups[last_j][-1][0]
                        j_y = self.groups[last_j][-1][1]
                        j_diagonal = j_y - j_x

                        if L_dict.peekitem()[1][0] < V[k] or abs(j_diagonal - diagonal) > 10:
                            L_dict[h_k] = (V[k], k)

                """
                if len(L) == 0 or L.max_item()[1][0] < V[k]:
                    L[h_k] = (V[k], k)
                """
                # L.insert(h_k, (V[k], k))

                j1_index, j1_key, j1_value = find_not_smaller_key(L_dict, h_k)

                if j1_index != None:
                    j1_index += 1
                    while j1_index < len(L_dict):
                        prev_j1_key = j1_key
                        prev_j1_value = j1_value
                        try:

                            j1_key = L_dict.iloc[j1_index]
                            j1_value = L_dict[j1_key]
                            if V[k] > j1_value[0]:
                                del L_dict[j1_key]
                            else:
                                j1_index += 1
                        except KeyError:
                            # print prev_j1_item
                            if V[k] > prev_j1_value[0]:
                                del L_dict[prev_j1_key]
                            break

        # print "DP finished"
        try:
            # print(back_track)
            # print(L_dict)
            current_j = V.index(max(V))
            score = V[current_j]
            # backtrack
            chain_index = []
            chain_index.append(current_j)

            while True:
                prev_j = back_track[current_j]
                if prev_j == -1:
                    break
                else:
                    current_j = prev_j
                    chain_index.append(current_j)

            optimal_chain = []
            length = 0
            for i in chain_index[::-1]:
                optimal_chain.extend(self.groups[i])
                length += sum([x[2] for x in self.groups[i]])

        except ValueError:
            optimal_chain = None
            length = 0

        return optimal_chain, length

    def rc_chain_groups(self):

        r = len(self.I)

        L_dict = SortedDict()
        V = [0] * len(self.groups)
        back_track = [-1] * len(self.groups)

        self.I.sort()

        for i in range(r):
            # I_list[i] is a start point, noticed we go through the chain from botton to top
            if self.I[i][2] == 0:
                k = self.I[i][1]

                l_k = self.groups[k][0][1] + self.groups[k][0][2]
                start_y = l_k
                start_x = self.I[i][0]
                end_x = self.groups[k][-1][0]
                end_y = self.groups[k][-1][1]
                diagonal = start_y + start_x
                # find largest h_j strictly smaller than l_k and also not off diagonal

                j_index, j_key, j_value = find_not_smaller_key(L_dict, l_k + 1)
                # print("j", j_value)
                while j_index != None:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = j_value[1]
                    j_x = self.groups[j][-1][0]
                    j_y = self.groups[j][-1][1]
                    d_x = start_x - j_x
                    d_y = j_y - start_y
                    # print("start",k, start_x, start_y)
                    # print(j, j_x, j_y)
                    # print(d_x, d_y)
                    if min(d_x, d_y) / max(d_x, d_y) > self.ratio:

                        prev_score = V[j]
                        break

                    else:
                        j_index, j_key, j_value = find_not_smaller_key(L_dict, j_y + 1)
                    # print(j)

                else:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = -1
                    prev_score = 0

                V[k] = prev_score + v_j
                back_track[k] = j

            # is a end point
            else:
                # print L
                k = self.I[i][1]
                h_k = self.groups[k][-1][1]

                l_k = self.groups[k][0][1] - self.groups[k][0][2]
                start_y = l_k
                start_x = self.I[i][0]
                diagonal = start_y + start_x

                j_index, j_key, j_value = find_not_bigger_key(L_dict, h_k)

                if j_index != None:
                    j = j_value[1]
                    V_j = j_value[0]

                    j_x = self.groups[j][-1][0]
                    j_y = self.groups[j][-1][1]
                    j_diagonal  = j_x + j_y

                    if V[k] > V_j or abs(j_diagonal - diagonal) > 10:
                        L_dict[h_k] = (V[k], k)

                else:

                    if len(L_dict) == 0:
                        L_dict[h_k] = (V[k], k)

                    else:
                        # print(L_dict.peekitem())
                        last_j = L_dict.peekitem(0)[1][1]
                        j_x = self.groups[last_j][-1][0]
                        j_y = self.groups[last_j][-1][1]
                        j_diagonal = j_y + j_x

                        if L_dict.peekitem(0)[1][0] < V[k] or abs(j_diagonal - diagonal) > 10:
                            L_dict[h_k] = (V[k], k)
                """
                if len(L) == 0 or L.max_item()[1][0] < V[k]:
                    L[h_k] = (V[k], k)
                """

                j1_index, j1_key, j1_value = find_not_bigger_key(L_dict, h_k)
                if j1_index != None:
                    j1_index -= 1
                    while j1_index >= 0:
                        prev_j1_key = j1_key
                        prev_j1_value = j1_value
                        try:
                            # print j1_index
                            j1_key = L_dict.iloc[j1_index]
                            j1_value = L_dict[j1_key]
                            if V[k] > j1_value[0]:
                                del L_dict[j1_key]

                            j1_index -= 1
                        except KeyError:
                            # print prev_j1_item
                            if V[k] > prev_j1_value[0]:
                                del L_dict[prev_j1_key]
                            break

        try:

            current_j = V.index(max(V))
            score =V[current_j]

            #current_j = max_value[1]
            # backtrack
            chain_index = []
            chain_index.append(current_j)

            while True:
                prev_j = back_track[current_j]
                if prev_j == -1:
                    break
                else:
                    current_j = prev_j
                    chain_index.append(current_j)

            optimal_chain = []
            length = 0
            for i in chain_index[::-1]:
                optimal_chain.extend(self.groups[i])
                length += sum([x[2] for x in self.groups[i]])

        except ValueError:
            optimal_chain = None
            length = 0

        return optimal_chain, length

    def chain_groups(self, accuracy=0.8, group_distance=None, rechain_threshold=5, span_coefficient=1.0, identity=400, release = 1.5):
        if group_distance == None:
            group_distance = statistical_bound_of_waiting_time(accuracy, 9)

        if len(self.groups) != 0:
            extend = None
            if self.forward:
                align, length = self.fw_chain_groups()

                if align:
                    # find left side extension length
                    if align[0][0] < align[0][1]:
                        left_extend = align[0][0] - align[0][2]
                        self.query_overlap_start = 0
                        self.target_overlap_start = align[0][1] - align[0][0] + align[0][2]
                    else:
                        left_extend = align[0][1] - align[0][2]
                        self.query_overlap_start = align[0][0] - align[0][1] + align[0][2]
                        self.target_overlap_start = 0

                    if self.query_len - align[-1][0] < self.target_len - align[-1][1]:
                        right_extend = self.query_len - align[-1][0]
                        self.query_overlap_end = self.query_len
                        self.target_overlap_end = align[-1][1] + self.query_len - align[-1][0]
                    else:
                        right_extend = self.target_len - align[-1][1]
                        self.query_overlap_end = align[-1][0] + self.target_len - align[-1][1]
                        self.target_overlap_start = self.target_len

                    # middle span
                    middle_extend = 0.5 * (abs(align[-1][0] - align[0][0] + align[0][2])) + 0.5 * abs(
                        align[-1][1] - align[0][1] + align[0][2])

                    extend = left_extend + middle_extend + right_extend

                    self.query_ali_start = align[0][0] - align[0][2]
                    self.query_ali_end = align[-1][0]
                    self.target_ali_start = align[0][1] - align[0][2]
                    self.target_ali_end = align[-1][1]

                """
                extend = 2 * group_distance + 0.5 * (
                    align[-1][0] - align[0][0] + align[0][2] + align[-1][1] - align[0][1] + align[0][2])
                """

            else:
                align, length = self.rc_chain_groups()

                if align:
                    # find left side extension length
                    if align[0][0] < self.target_len - align[0][1]:
                        left_extend = align[0][0] - align[0][2]
                        self.query_overlap_start = 0
                        self.target_overlap_end = align[0][1] + align[0][0]
                    else:
                        left_extend = self.target_len - align[0][1]
                        self.query_overlap_start = align[0][0] - self.target_len + align[0][1]
                        self.target_overlap_end = self.target_len

                    if self.query_len - align[-1][0] < align[-1][1]:
                        right_extend = self.query_len - align[-1][0]
                        self.query_overlap_end = self.query_len
                        self.target_overlap_start = align[-1][1] - self.query_len + align[-1][0]
                    else:
                        right_extend = align[-1][1]
                        self.query_overlap_end = align[-1][0] + align[-1][1]
                        self.target_overlap_start = 0

                    # middle span
                    middle_extend = 0.5 * (abs(align[-1][0] - align[0][0] + align[0][2])) + 0.5 * abs(
                        align[0][1] - align[-1][1] + align[0][2])

                    extend = left_extend + middle_extend + right_extend
                    self.query_ali_start = align[0][0] - align[0][2]
                    self.query_ali_end = align[-1][0]
                    self.target_ali_start = align[-1][1]
                    self.target_ali_end = align[0][1] + align[0][2]

                """
                extend = 2 * group_distance + 0.5 * (
                    align[-1][0] - align[0][0] + align[0][2] + align[0][1] - align[-1][1] + align[0][2])
                """

            # compare align region of two sequnces to make sure they are close to clollinear
            query_ali_len = abs(self.query_ali_end - self.query_ali_start)
            target_ali_len = abs(self.target_ali_end - self.target_ali_start)
            ratio = min(query_ali_len, target_ali_len) / max(query_ali_len, target_ali_len)

            """
            print()
            print(align)
            print("number of groups", len(align))
            print("length", length)
            print(length / 9)
            print(extend / float(span_coefficient * group_distance))
            print(ratio)
            """

            if ratio >= self.ratio and extend:
                if (length > identity and extend / float(span_coefficient * group_distance) <= float(release * length) / 9) or (
                            extend / float(span_coefficient * group_distance) <= float(
                            length) / 9 and length > rechain_threshold * 9):
                    self.chain_align = align
                    self.aligned = True
                    self.aligned_base = length
                    self.overlap_length = extend

    def set_ratio(self, ratio):
        self.ratio = ratio

    def plot_hits(self):
        print("groups", self.groups)
        print(self.query, self.target)
        plt.figure()

        if self.aligned:
            print("chains", self.chain_align)
            coordinates = list(map(list, zip(*self.chain_align)))
            plt.scatter(coordinates[0], coordinates[1], s=40, edgecolors="black", linewidths=5)

        for group in self.groups:
            # print chain
            chain_coor = list(map(list, zip(*group)))
            # print(chain_coor)
            plt.scatter(chain_coor[0], chain_coor[1], s=40)

        plt.show()


if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("k", help="kmer", type=int)
    parser.add_argument("accuracy", help="accuracy of the reads", type=float)
    # parser.add_argument("gap", help="gap rate", type=float)
    parser.add_argument("rechain", help="length of kmer", type=int)
    parser.add_argument("span_coefficient", help="must be non-zero float", type=float)
    parser.add_argument("identical_base",
                        help="if the number of identical base meet this threshold the output will always be reported",
                        type=int)
    parser.add_argument("--ratio", help="ratio of two overlap region threshold", type=float, default=0.5)
    parser.add_argument("--size", help="group size threshold", type=int, default=12)
    parser.add_argument("--large", help="release threhold for large overlap", type=float, default=1.5)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    file = args.input
    L = statistical_bound_of_waiting_time(args.accuracy, args.k)

    with open(file) as f:
        lines = f.readlines()

    for line in lines:
        group_hit = GroupHit(line, args.size)
        group_hit.set_ratio(args.ratio)

        # print group_hit.groups

        group_hit.chain_groups(accuracy=args.accuracy, group_distance=L, rechain_threshold=args.rechain,
                               span_coefficient=args.span_coefficient, identity=args.identical_base, release= args.large)
        # print group_hit.chain_align
        if group_hit.aligned:
            output_str = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                group_hit.query, group_hit.target, str(group_hit.aligned), group_hit.aligned_base,
                group_hit.overlap_length, group_hit.query_ali_start, group_hit.query_ali_end,
                group_hit.target_ali_start, group_hit.target_ali_end, group_hit.query_overlap_start,
                group_hit.query_overlap_end, group_hit.target_overlap_start, group_hit.target_overlap_end
            )

            print(output_str)

        # group_hit.plot_hits()