class Duplication:
    def __init__(self, pos_in_ref, pos_in_qry, unit_len, n_dup, is_reverse):
        self.pos_in_ref = pos_in_ref
        self.pos_in_qry = pos_in_qry
        self.unit_len = unit_len
        self.n_dup = n_dup
        self.is_reverse = is_reverse

    def print_dup(self):
        print(f"pos_in_ref={self.pos_in_ref}    "
              f"pos_in_qry={self.pos_in_qry}    "
              f"unit_len={self.unit_len}    "
              f"n_dup={self.n_dup}  "
              f"is_reverse={self.is_reverse}\n")

    def dec(self):
        self.n_dup -= 1


def get_reverse_complement_of(s):
    complement_of = {'C': 'G', 'A': 'T', 'G': 'C', 'T': 'A'}
    s_rc = ''
    for ch in reversed(s):
        s_rc += complement_of[ch]
    return s_rc

def get_first_different(s1, s2):
    i = 0
    while s1[i] == s2[i]:
        i += 1
    return i


class DNAMatcher:

    def __init__(self, query, reference):
        self.qry = query
        self.ref = reference

    @classmethod
    def read_file(cls, qry_path, ref_path):
        qry = ''
        ref = ''
        with open(qry_path, 'r') as qry_file:
            qry = qry_file.read().replace('\n', '')
        with open(ref_path, 'r') as ref_file:
            ref = ref_file.read().replace('\n', '')

        return cls(qry, ref)

    def get_pos_in_ref(self, l):
        pos_in_ref = {}
        l_ref = len(self.ref)

        for start in range(0, l_ref - l + 1):
            cur_str = self.qry[start: start + l]
            cur_str_rc = get_reverse_complement_of(cur_str)

            if cur_str not in pos_in_ref:
                pos_in_ref[cur_str] = start + l
            if cur_str_rc not in pos_in_ref:
                pos_in_ref[cur_str_rc] = start + l

        return pos_in_ref

    def get_n_duplicated(self, target_str, start, l):
        n_duplicated = 0

        while True:
            if start + l >= len(self.qry):
                break
            if target_str != self.qry[start: start + l]:
                break
            start = start + l
            n_duplicated += 1

        return n_duplicated + 1


    def match(self):
        """
        思路：
        对于每一个长度l(l=1,...,l_qry)，我们：
        1. 记录所有长度为l的子串在ref中出现的位置，记录为pos_in_ref(参照文档中给出的格式，位置用该子串的结束位置，即串联重复的开始位置表示)
           同时，每个子串的反向互补串也需要记录。
        2. 扫描qry中所有长度为l的子串s。如果该子串(或其反向互补串)在pos_in_ref中有记录，那么计算s在qry中的重复次数

        需要注意的一点是，由于连续重复的周期性，重复序列可能会有多种不同的理解，导致
        一次真重复会被检测出很多个伪重复
        为了避免这一点，我们优先寻找最长的重复序列，并且为已经被长重复序列覆盖到的区域打上标记

        """

        len_qry = len(self.qry)

        duplications = []
        has_probed = [False] * len_qry

        first_diff = get_first_different(self.qry, self.ref)

        for l in range(200, 20, -1):
            pos_in_ref = self.get_pos_in_ref(l)
            start = 0
            while start <= len_qry - l:

                if has_probed[start]:
                    start += 1
                    continue

                cur_str = self.qry[start: start + l]
                if cur_str not in pos_in_ref:
                    start += 1
                    continue

                cur_pos_in_ref = pos_in_ref[cur_str]

                if cur_pos_in_ref > first_diff:
                    start += 1
                    continue

                is_reverse = (cur_str != self.ref[cur_pos_in_ref - l: cur_pos_in_ref])

                n_dup = self.get_n_duplicated(cur_str, start + l, l)
                if n_dup == 1:
                    start += 1
                    continue
                duplications.append(Duplication(cur_pos_in_ref, start, l, n_dup, is_reverse))
                for j in range(start, min(start + l * n_dup, len_qry)):
                    has_probed[j] = True
                start = start + l * n_dup
        duplications[len(duplications) - 1].dec()
        return duplications


matcher = DNAMatcher.read_file("query.txt", "reference.txt")

dups = matcher.match()
for dup in dups:
    dup.print_dup()
