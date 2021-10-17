from math import exp
import numpy as np
import matplotlib.pyplot as plt

########################################################################
# Constants
########################################################################

# [V]
Vthn = 2.00E-01
# [V]
Vthp = 2.32E-01
# [A/K^2]
a = 2.02E-07
# [K/V]
b = 9.21E+03
# [K/V]
c = 8.90E+02
#
g = 2.00E-01
# [Farad]
C = 2.50E-01
#
s = 5.00E-04
#
Omega = 5.00E-01
# [V]
Vmin = 1.00E+00
# [MHz]
fmin = 8.00E+00
# [MHz]
fmax = 2.00E+02
# [MHz/V]
gamma = 2.00E+03
# [s] - Actuation Period
p = 1.00E+00
# [s] - Time to switch frequency
Tfreq = 1.00E-06
# [s] - Time from active to clock gate
Tact_clock = 1.00E-06
# [s] - Time from clock gate to active
Tclock_act = 1.00E-06
# [s] - Time from active to power gate
Tact_power = 1.00E-05
# [s] - Time from active to power gate
Tpower_act = 1.00E-03
# [Ah] - Battery
Bat = 5.40E+00
# [J] Energy available
Eavailable = 1.5 * Bat * 3600

def vdd(f):
    return Vmin + (f - fmin) / gamma

def calc_days_of_life(e_consumed_per_period):
    return Eavailable / e_consumed_per_period / 60 / 60 / 24

def calculate_total_energy(
        # Active frequency
        f_active,
        # Environment temperature
        temperature,
        # False = constant frequency
        # True = change to fmin before sleeping
        low_freq_sleep=True,
        # None
        # 'clock'
        # 'power'
        sleep_strategy='clock',
        # Print reports
        print_report=False,
        # Tactive = Tactive * active_cycle_optimization
        # Example: optmization = 0.2, Tactive = 1s -> new Tactive = 0.2s
        active_cycle_optimization=1,
        debug=False):
    
    T = temperature
    f = f_active

    ########################################################################
    #Auxiliary calculations =  #
    ########################################################################

    # - First exponential of Ps formula
    Exp1 = exp(-b * Vthp / T)
    # - Second exponential of Ps formula
    Exp2 = exp(-b * Vthn / T)

    ########################################################################
    # Formulas
    ########################################################################

    if sleep_strategy == None:
        Tactive = p
        Tto_active = 0
        Tto_sleep = 0
        Tsleep = 0
    else:
        # [s] - Time processor is active = p * 1 / (f / fmin)
        Tactive = p * 1 / (f / fmin)
        Tactive = Tactive * active_cycle_optimization
        # [s] - "Delay till the processor is active.
        # Selected between: Tclock_act, Tpower_act"
        Tto_active = Tclock_act if sleep_strategy == 'clock' else Tpower_act
        # [s] - "Delay till the processor is asleep.
        # Selected between: Tact_clock Tact_power"
        Tto_sleep = Tact_clock if sleep_strategy == 'clock' else Tact_power
        # [s] - "Time processor is in sleep mode.
        # If frequency does not change:
        #   Tsleep = p - (Tactive + Tto_active + Tto_sleep)"
        # If frequency changes: 
        #   Tsleep = p - (Tactive + Tto_active + Tto_sleep + 2 * Tfreq)"
        if low_freq_sleep:
            Tsleep = max(0, 
                p - (Tactive + Tto_active + Tto_sleep + 2 * Tfreq))
        else:
            Tsleep = max(0, p - (Tactive + Tto_active + Tto_sleep))

    # [W]
    Ps_fmin = a * vdd(fmin) * T * T * (Exp1 + Exp2) * exp(c * vdd(fmin) / T)
    # [W]
    Ps_fact = a * vdd(f) * T * T * (Exp1 + Exp2) * exp(c * vdd(f) / T)
    # [W]
    Pa = C * (vdd(f) ** 2) * f * s + Omega * Ps_fact

    # [J] 
    # Total energy consumed during transition from active to sleep + sleep to active"
    Etransition = Pa * (Tto_active + Tto_sleep)
    # [J] - Total energy consumed while active
    if low_freq_sleep:
        Echange_freq = Pa * 2 * Tfreq
    else:
        Echange_freq = 0
    # [J] - Total energy consumed while active
    Eactive = Pa * Tactive
    # [J] - Total energy consumed while asleep
    if low_freq_sleep:
        if sleep_strategy == 'clock':
            Esleep = Ps_fmin * Tsleep
        else:
            Esleep = g * Ps_fmin * Tsleep
    else:
        if sleep_strategy == 'clock':
            Esleep = Ps_fact * Tsleep
        else:
            Esleep = g * Ps_fact * Tsleep

    # [J] - Total energy consumed during 1 period
    Etotal = Eactive + Etransition + Esleep + Echange_freq

    days_active = calc_days_of_life(Etotal)

    if sleep_strategy == None:
        sleep_text = 'Always on'
    else:
        sleep_text = sleep_strategy + 'gated'

    if low_freq_sleep:
        freq_policy = 'Sleep @ 8MHz'
    else:
        freq_policy = 'No frequency change'

    if print_report:
        print(f'Energy report for:')
        print(f"f active:         {f} MHz")
        print(f"Temperature:      {T - 273}" + u'\N{DEGREE SIGN}' + 'C')
        print(f"Sleep strategy:   {sleep_text}")
        print(f"Frequency policy: {freq_policy}")
        print('-' * 40)
        print(f'Period time profiling (percentage to total period time):')
        print(f"Tsleep:         {Tsleep:.3E} s - {100 * Tsleep/p:.2f}%")
        print(f"Tto_sleep:      {Tto_sleep:.3E} s - {100 * Tto_sleep/p:.2f}%")
        print(f"Tto_active:     {Tto_active:.3E} s - {100 * Tto_active/p:.2f}%")
        print(f"Tactive:        {Tactive:.3E} s - {100 * Tactive/p:.2f}%")
        if low_freq_sleep:
            print(f"Tchange_freq:        {Tfreq:.3E} s - {100 * Tfreq/p:.2f}%")
        if debug:
            print('-' * 40)
            print(f'Power dissipation:')
            print(f"Ps_fmin:        {Ps_fmin:.3E} W")
            print(f"Ps_fact:        {Ps_fact:.3E} W")
            print(f"Pa:             {Pa:.3E} W")
        print('-' * 40)
        print(f'Energy Consumption:')
        print(f"E transition:   {1000 * Etransition:.3E} mJ - {100 * Etransition/Etotal:.2f}%")
        print(f"E change freq:  {1000 * Echange_freq:.3E} mJ - {100 * Echange_freq/Etotal:.2f}%")
        print(f"E sleep:        {1000 * Esleep:.3E} mJ - {100 * Esleep/Etotal:.2f}%")
        print(f"E active:       {1000 * Eactive:.3E} mJ - {100 * Eactive/Etotal:.2f}%")
        print(f"E total:        {1000 * Etotal:.3E} mJ")
        print(f"\n >> Days of life: {days_active:.2f}")
        print('\n' + '-/' * 20 + '-')

    data = {
            'Tsleep': Tsleep,
            'Tto_sleep': Tto_sleep,
            'Tto_active': Tto_active,
            'Tactive': Tactive,
            'Etransition': Etransition,
            'Echange_freq': Echange_freq,
            'Esleep': Esleep,
            'Eactive': Eactive,
            'Etotal': Etotal
           }
    return Etotal

def annot_min(x,y, t_index, total_t, ax=None):
    text= f"{y:.3f}mJ @ {x}MHz"
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="data",
              arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(x, y), xytext=(x + 190 - 18 * t_index, 1.1 * y), **kw)

def generate_plots(
        title='Energy Consumption by period',
        filename='',
        strategy='power',
        low_freq=True,
        annot_min_point=False,
        active_cycle_optimization=1,
        freqs_range=(8, 201, 1),
        temps_range=(273 - 25, 273 + 66, 15)):

    freqs = np.arange(*freqs_range)
    temps = np.arange(*temps_range)
    
    best_results = []

    plt.subplot(1, 2, 1)

    for tindex, t in enumerate(temps):
        data = [calculate_total_energy(f, t, 
                    sleep_strategy=strategy,
                    active_cycle_optimization=active_cycle_optimization,
                    low_freq_sleep=low_freq)
                for f in freqs]
        best_consumption = min(data)
        best_freq = np.argmin(data) + fmin
        # Converting to mJ for better plot
        data = [x * 1000 for x in data]
        best_results.append({'temperature': t,
                             'freq': best_freq,
                             'consumption': best_consumption})
        plt.plot(freqs, data, label=f'T={t-273}' + u'\N{DEGREE SIGN}' + 'C')
        if annot_min_point:
            annot_min(best_freq, 1000 * best_consumption, tindex, len(temps))
    
    plt.legend()
    plt.title(title)
    plt.xlabel('Active frequency (MHz)')
    plt.ylabel('Energy Consumption (mJ)')

    # Dayd of life plot
    plt.subplot(1, 2, 2)
    days_of_live = [calc_days_of_life(x['consumption'])
                    for x in best_results]
    x_axis = [('{}' + u'\N{DEGREE SIGN}' + 'C\n' + '{}MHz')
                    .format(x['temperature'] - 273, x['freq'])
                for x in best_results]
    plt.bar(x_axis, days_of_live)
    plt.title('Days of Life @ Best Frequency')
    plt.xlabel('Temperature (' + u'\N{DEGREE SIGN}' + 'C)')
    plt.ylabel('Days of life')
    

    plt.tight_layout()
    figure = plt.gcf()
    figure.set_size_inches(15, 12)
    if filename:
        plt.savefig(filename + '.png', dpi=100)
    plt.show()

    return best_results

########################################################################
# Baseline scenario - Constant frequency, no sleep

# calculate_total_energy(8, 273 - 20, sleep_strategy=None, print_report=True)
# calculate_total_energy(8, 273 + 20, sleep_strategy=None, print_report=True)
# calculate_total_energy(8, 273 + 65, sleep_strategy=None, print_report=True)
# generate_plots(title='Consumption by Duty Cycle',
#         filename='always_on',
#         strategy=None)

########################################################################
# Clock gated sleep optimization

# e_no_low = calculate_total_energy(8, 273 + 25,
#         sleep_strategy='clock', low_freq_sleep=False,
#         print_report=True)
# e_with_low = calculate_total_energy(8, 273 + 25,
#         sleep_strategy='clock', print_report=True)
# generate_plots(title='Consumption by Duty Cycle',
#         low_freq=False,
#         filename='clock_gated_const_freq',
#         annot_min_point=True,
#         strategy='clock',
# generate_plots(title='Consumption by Duty Cycle',
#         low_freq=True,
#         filename='clock_gated_low_freq',
#         annot_min_point=True,
#         strategy='clock')

########################################################################
# Power gated sleep optimization

# e_no_low = calculate_total_energy(25, 273 - 25,
#         sleep_strategy='power',
#         low_freq_sleep=False,
#         print_report=True)
# e_with_low = calculate_total_energy(37, 273 + 20,
#         low_freq_sleep=False,
#         sleep_strategy='power', print_report=True)
# e_with_low = calculate_total_energy(49, 273 + 65,
#         low_freq_sleep=False,
#         sleep_strategy='power', print_report=True)
# generate_plots(title='Consumption by Duty Cycle',
#         low_freq=False,
#         filename='power_gated_const_freq',
#         annot_min_point=True,
#         strategy='power')

# e_no_low = calculate_total_energy(27, 273 - 25,
#         sleep_strategy='power',
#         print_report=True)
# e_with_low = calculate_total_energy(43, 273 + 20,
#         sleep_strategy='power', print_report=True)
# e_with_low = calculate_total_energy(62, 273 + 65,
#         sleep_strategy='power', print_report=True)
# generate_plots(title='Consumption by Duty Cycle',
#         low_freq=True,
#         filename='power_gated_low_freq',
#         annot_min_point=True,
#         strategy='power')

########################################################################
# Power gated sleep + DMA optimization

# generate_plots(title='Consumption by Duty Cycle',
#         low_freq=True,
#         filename='dma_20',
#         active_cycle_optimization=0.8,
#         annot_min_point=True,
#         strategy='power')
e_no_low = calculate_total_energy(27, 273 - 25,
        sleep_strategy='power',
        active_cycle_optimization=0.8,
        print_report=True)
e_with_low = calculate_total_energy(43, 273 + 20,
        active_cycle_optimization=0.8,
        sleep_strategy='power', print_report=True)
e_with_low = calculate_total_energy(62, 273 + 65,
        active_cycle_optimization=0.8,
        sleep_strategy='power', print_report=True)
