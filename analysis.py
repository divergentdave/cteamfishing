#!/usr/bin/env python3
import collections
import fractions
import functools
import operator


class FishingGameAnalysis:
    FULL_NET = (12, 10, 10, 8, 6, 4)

    def __init__(self):
        self.optimal_choices = {}
        self.net_pdf_memo = {}
        self.subset_sum_memo = {}

    def net_pdf(self, net):
        """
        The input of this function is a tuple of which dice are in the net
        before rolling.
        The output of this function is a dictionary representing the
        probability distribution of playing out the game from this roll until
        the end of the cast. The probability of getting some number of fish
        over the rest of the cast can be obtained by looking up that number
        as a key in the dictionary. (The values in the dictionary should thus
        have a sum of 1)
        """
        if net in self.net_pdf_memo:
            return self.net_pdf_memo[net]

        pdf_accum = collections.defaultdict(lambda: fractions.Fraction(0))
        if len(net) == 0:
            pdf_accum[0] = fractions.Fraction(1)
            self.net_pdf_memo[net] = pdf_accum
            return pdf_accum
        denominator = 20 * self.roll_combinations(net)
        one_over_denom = fractions.Fraction(1, denominator)

        def catch_update_pdf(next_net, factor):
            temp_pdf = self.net_pdf(next_net)
            for fish, probability in temp_pdf.items():
                pdf_accum[fish + 1] += probability / factor

        # special case: 1
        for dice in self.roll_gen(net):
            if 1 in dice:
                next_nets = list(self.nat1_next_nets(net, dice))
                next_net = self.choose_next_net(next_nets)
                catch_update_pdf(next_net, denominator)
            else:
                pdf_accum[0] += one_over_denom

        # normal case: 2-19
        for target in range(2, 20):
            for dice in self.roll_gen(net):
                next_nets = list(self.add_sub_next_nets(net, dice, target))
                if next_nets:
                    next_net = self.choose_next_net(next_nets)
                    catch_update_pdf(next_net, denominator)
                else:
                    pdf_accum[0] += one_over_denom

        # special case: 20
        next_net = self.choose_next_net(list(self.nat20_next_nets(net)))
        catch_update_pdf(next_net, 20)

        assert sum(pdf_accum.values()) == 1
        self.net_pdf_memo[net] = pdf_accum
        return pdf_accum

    def roll_combinations(self, net):
        """
        Returns how many permutations can be rolled with the dice in the net.
        For simplicity, it's assumed the two D10s are distinguished from each
        other.
        """
        return functools.reduce(operator.mul, net)

    def roll_gen(self, net):
        """
        Given a tuple of dice in the net, yields each possible roll as a tuple.
        """
        if len(net) == 1:
            for i in range(1, net[0] + 1):
                yield (i,)
        else:
            for partial in self.roll_gen(net[:-1]):
                for i in range(1, net[-1] + 1):
                    yield partial + (i,)

    def nat1_next_nets(self, net, dice):
        """
        Given the current net, and the results of their roll, yields the next
        possible nets to catch a fish with a target of 1.
        """
        dedup = set()
        for i in range(len(net)):
            if dice[i] == 1:
                next_net = net[:i] + net[i + 1:]
                if next_net not in dedup:
                    dedup.add(next_net)
                    yield next_net

    def nat20_next_nets(self, net):
        "Given the current net, yields the next possible nets on a natural 20."
        dedup = set()
        for i in range(len(net)):
            next_net = net[:i] + net[i + 1:]
            if next_net not in dedup:
                dedup.add(next_net)
                yield next_net

    def subset_sums_cached(self, numbers, target):
        """
        Given a set of numbers and a target number, returns a sequence of
        tuples of indices into the input numbers, so that for any such tuple,
        sum(numbers[i] for i in tuple) will be equal to the target number.
        Internally, the results are cached, and the numbers are sorted to take
        advantage of symmetries when caching results, but this doesn't affect
        the indices that are returned.
        """
        numbers_and_indices = zip(numbers, range(len(numbers)))
        numbers_and_indices_sorted = sorted(numbers_and_indices)
        numbers_sorted = tuple(tup[0] for tup in numbers_and_indices_sorted)

        def subset_sums(remaining_metaindices):
            if remaining_metaindices:
                for partial_sum, partial_subset in subset_sums(
                        remaining_metaindices[1:]):
                    number = numbers_sorted[remaining_metaindices[0]]
                    yield (
                        partial_sum - number,
                        (remaining_metaindices[0],) + partial_subset
                    )
                    yield partial_sum, partial_subset
                    yield (
                        partial_sum + number,
                        (remaining_metaindices[0],) + partial_subset
                    )
            else:
                yield 0, ()

        if numbers_sorted not in self.subset_sum_memo:
            self.subset_sum_memo[numbers_sorted] = {
                n: [] for n in range(2, 20)
            }
            for total, metaindices in subset_sums(list(range(len(numbers)))):
                if total >= 2 and total <= 19:
                    self.subset_sum_memo[numbers_sorted][total].append(
                        metaindices
                    )
        for metaindices in self.subset_sum_memo[numbers_sorted][target]:
            yield [
                numbers_and_indices_sorted[metaindex][1]
                for metaindex in metaindices
            ]

    def add_sub_next_nets(self, net, dice, target):
        dedup = set()
        for indices in self.subset_sums_cached(dice, target):
            next_net = tuple(
                net[i]
                for i in range(len(net))
                if i not in indices
            )
            if next_net not in dedup:
                dedup.add(next_net)
                yield next_net

    def choose_next_net(self, next_net_choices):
        if len(next_net_choices) == 1:
            return next_net_choices[0]
        normalized = tuple(sorted(next_net_choices))
        if normalized in self.optimal_choices:
            return self.optimal_choices[normalized]
        evs = [self.pdf_to_ev(self.net_pdf(net)) for net in normalized]
        optimal_index = max(list(range(len(evs))), key=lambda i: evs[i])
        optimal_net = normalized[optimal_index]
        self.optimal_choices[normalized] = optimal_net
        return optimal_net

    def pdf_to_ev(self, pdf):
        return sum(fish * probability for fish, probability in pdf.items())


def main():
    global analysis
    analysis = FishingGameAnalysis()

    def flatten_pdf(pdf_dict):
        return [float(pdf_dict[i]) for i in range(len(analysis.FULL_NET) + 1)]

    print(flatten_pdf(analysis.net_pdf(analysis.FULL_NET)))


if __name__ == "__main__":
    main()
