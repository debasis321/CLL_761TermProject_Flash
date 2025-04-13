from graphviz import Digraph

# Create a colorful flowchart including corrected phase check position
dot = Digraph(comment='Peng-Robinson Flash Algorithm Corrected', format='png')
dot.attr(rankdir='TB', size='5,5')
dot.attr(ranksep='0.2')    # Reduce vertical spacing
dot.attr(dpi='2400')
concentrate='true'

# Node styles
style_common = {'style': 'filled', 'fontname': 'Helvetica', 'fontsize': '10'}
colors = {
    'start': 'lightblue',
    'data': 'violet',
    'estimate': 'gold',
    'calc': 'orange',
    'solve': 'lightgreen',
    'update': 'lightcoral',
    'decision': 'lightgray',
    'end': 'lightskyblue',
    'state': 'lightyellow',
    'adjuster': 'lightpink'
}

# Add nodes with updated flow
dot.node('A', 'Start: Feed I/P (T, P, z_i)', fillcolor=colors['start'], **style_common)
dot.node('A1', 'Load Thermodynamic Data\n(T_c, P_c, ω, MW, k_ij, Cp(T))', fillcolor=colors['data'], **style_common)
dot.node('A2', 'Calculation of Residual Property\n(T_r, P_r)', fillcolor=colors['calc'], **style_common)
dot.node('B', 'Estimate K_i (Wilson)', fillcolor=colors['estimate'], **style_common)
dot.node('C', 'Estimate β, x_i, y_i', fillcolor=colors['estimate'], **style_common)
dot.node('D', 'Calc EOS params: a_i, b_i, A, B', fillcolor=colors['calc'], **style_common)
dot.node('E', 'Solve \n(Cubic Euation for Z_l, Z_v)', fillcolor=colors['solve'], **style_common)
dot.node('E2', 'Check Phase Type\n(Liquid / Vapor / Two-Phase)\n(K_i based)', shape='diamond', fillcolor=colors['decision'], **style_common)
dot.node('F', 'Compute φ_i^L, φ_i^V', fillcolor=colors['calc'], **style_common)
dot.node('G', 'Update K_i = φ_i^L / φ_i^V', fillcolor=colors['update'], **style_common)
dot.node('H', 'Solver \n(Rachford-Rice for β)', fillcolor=colors['solve'], **style_common)
dot.node('H2', 'Check Phase Type\n(Liquid / Vapor / Two-Phase)\n(K_i based)', shape='diamond', fillcolor=colors['decision'], **style_common)
dot.node('I', 'Check Convergence\n(x_i, y_i, beta)', shape='diamond', fillcolor=colors['decision'], **style_common)
dot.node('J', 'Return x_i, y_i, β', fillcolor=colors['end'], **style_common)
dot.node('K1', 'Single Phase: Liquid', fillcolor=colors['end'], **style_common)
dot.node('K2', 'Single Phase: Vapor', fillcolor=colors['end'], **style_common)

dot.node('State1', 'Stream I/P State \n (T, P)', shape='box', fillcolor=colors['data'], **style_common)
dot.node('Solv1', 'Solver\n (Heat Balance)', fillcolor=colors['solve'], **style_common)
dot.node('Target', 'Target I/P\n(dP, dT, dQ, beta) ', shape='box', fillcolor=colors['adjuster'], **style_common)

dot.node('State2', 'Stream O/P State \n (T, P, x_i, y_i)', shape='box', fillcolor=colors['data'], **style_common)
dot.node('CalcRH', 'Calculate Resid. H ', fillcolor=colors['calc'], **style_common)
dot.node('IdealH', 'Ideal H', fillcolor=colors['calc'], **style_common)
dot.node('ActualH', 'Actual H', fillcolor=colors['calc'], **style_common)
dot.node('O/P', 'End: Q_Flash O/P', fillcolor=colors['start'], **style_common)

dot.edge('Target', 'Solv1', )
dot.edge('Solv1', 'State1', )
dot.edge('State1', 'A2', label='Trigger Flash Calculation')
dot.edge('A', 'State1', label='Input', linetype='dashed')
# Edges with corrected flow
dot.edge('A', 'A1')
dot.edge('A1', 'A2')
dot.edge('A2', 'B')
dot.edge('B', 'C')
dot.edge('C', 'D')
dot.edge('D', 'E')
dot.edge('E', 'E2')
dot.edge('E2', 'F', label='Two-Phase')
dot.edge('F', 'G')
dot.edge('G', 'H')
dot.edge('H', 'H2')
dot.edge('H2', 'K1', label='Liquid')
dot.edge('H2', 'K2', label='Vapor')
dot.edge('E2', 'K1', label='Liquid')
dot.edge('E2', 'K2', label='Vapor')
dot.edge('H2', 'I', label='Two-Phase')
dot.edge('I', 'C', label='No')
dot.edge('I', 'J', label='Yes')
dot.edge('J', 'State2')
dot.edge('K1', 'State2')
dot.edge('K2', 'State2')
dot.edge('State2', 'CalcRH')
dot.edge('A1', 'IdealH', label='Cp Coefficients')
dot.edge('State2', 'IdealH', label='T')
dot.edge('IdealH', 'ActualH')
dot.edge('CalcRH', 'ActualH')
dot.edge('ActualH', 'Solv1', label='H')
dot.edge('Solv1', 'O/P')

# Output path
corrected_flowchart_path = "flowChart.png"
dot.render(corrected_flowchart_path[:-4], format='png', cleanup=True)

corrected_flowchart_path
