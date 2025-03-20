import matplotlib.pyplot as plt
import numpy as np
import subprocess
import re
import os
from matplotlib.ticker import MaxNLocator
import pandas as pd
import seaborn as sns

class LDPCVisualizer:
    def __init__(self):
        # Parameters from C code
        self.num_monte_carlo = 0
        self.max_iterations = 0
        self.matrix_file = ""
        self.snr_value = 0.0
        self.n_vc = 0
        self.n_cv = 0
        self.offset_factor = 0.0
        self.max_operations = 0
        
        # Simulation results
        self.iterations_data = []
        self.fer_data = []
        self.ber_data = []
        self.syndrome_data = []
        self.frame_count_data = []
        
        # Decoding process data (for visualization)
        self.iter_llr_values = []  # Store LLR values at each iteration
        self.v_to_c_messages = []  # Variable to check node messages
        self.c_to_v_messages = []  # Check to variable node messages
        
    def parse_c_output(self, output_text):
        """Parse the output from the C simulation to extract relevant data."""
        # Extract simulation parameters
        params_match = re.search(r'NbMonteCarlo\s+:\s+(\d+).*?NbIterMax\s+:\s+(\d+).*?FileMatrix\s+:\s+([\w\.\/]+).*?snrValue\s+:\s+([\d\.]+).*?n_vc\s+:\s+(\d+).*?n_cv\s+:\s+(\d+).*?Offset\s+:\s+([\d\.]+).*?NbOper\s+:\s+(\d+)', 
                               output_text, re.DOTALL)
        
        if params_match:
            self.num_monte_carlo = int(params_match.group(1))
            self.max_iterations = int(params_match.group(2))
            self.matrix_file = params_match.group(3)
            self.snr_value = float(params_match.group(4))
            self.n_vc = int(params_match.group(5))
            self.n_cv = int(params_match.group(6))
            self.offset_factor = float(params_match.group(7))
            self.max_operations = int(params_match.group(8))
            
        # Extract frame results from progress reports
        frame_results = re.findall(r'<\d+>\s+FER=\s+(\d+)\s+/\s+(\d+)\s+=\s+([\d\.]+)\s+BER=\s+(\d+)\s+/\s+x\s+=\s+([\d\.]+)\s+avr_it=([\d\.]+)', output_text)
        
        for result in frame_results:
            frame_errors, frame_count, fer, bit_errors, ber, avg_iter = result
            self.fer_data.append(float(fer))
            self.ber_data.append(float(ber))
            self.iterations_data.append(float(avg_iter))
            self.frame_count_data.append(int(frame_count))
        
    def run_simulation(self, executable_path, args):
        """Run the C simulation with the given arguments and capture the output."""
        cmd = [executable_path] + args
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"Error running simulation: {stderr}")
            return False
        
        self.parse_c_output(stdout)
        return True
        
    def generate_synthetic_data(self):
        """Generate synthetic data for visualization when not running actual simulation."""
        # Code parameters
        n = 12  # Code length
        m = 6   # Number of check nodes
        k = n - m  # Information length
        gf = 64  # Galois field size
        
        # Generate synthetic codeword and noisy data
        codeword = np.random.randint(0, 2, size=(n, 6))  # Binary representation of GF symbols
        received = codeword.copy()
        
        # Add some noise (bit flips)
        noise_positions = np.random.choice(n*6, size=int(n*6*0.1), replace=False)
        for pos in noise_positions:
            i, j = pos // 6, pos % 6
            received[i, j] = 1 - received[i, j]
            
        # Create synthetic Tanner graph
        # Variable nodes connect to check nodes based on parity check matrix
        v_nodes = n
        c_nodes = m
        connections = []
        
        # Create random connections for the Tanner graph
        for c in range(c_nodes):
            # Each check node connects to some variable nodes
            connected_v = np.random.choice(v_nodes, size=3, replace=False)
            for v in connected_v:
                connections.append((v, c))
                
        # Synthetic LLR values for iterations
        iterations = 5
        llr_data = []
        for iter in range(iterations):
            # Generate LLR values that converge over iterations
            # For each symbol, generate LLR values for all GF elements
            iter_llr = np.zeros((n, gf))
            
            for i in range(n):
                # The true symbol has the most negative LLR (more likely to be 1)
                true_symbol = np.random.randint(0, gf)
                
                # Generate random LLRs for all symbols, with bias toward the true symbol
                # LLRs become more discriminative (larger magnitude) with each iteration
                base_llrs = np.random.normal(5.0, 2.0/(iter+1), size=gf)
                
                # Make the true symbol have strongly negative LLR
                base_llrs[true_symbol] = -10.0 * (iter + 1)
                
                # Normalize LLRs 
                base_llrs = base_llrs - np.min(base_llrs)
                iter_llr[i, :] = base_llrs
                
            llr_data.append(iter_llr)
            
        return {
            'codeword': codeword,
            'received': received,
            'v_nodes': v_nodes,
            'c_nodes': c_nodes,
            'connections': connections,
            'llr_data': llr_data,
            'iterations': iterations
        }
        
    def plot_tanner_graph(self, ax, v_nodes, c_nodes, connections):
        """Plot the Tanner graph representation of the LDPC code."""
        # Position variable nodes at the top
        v_pos_x = np.linspace(0, 1, v_nodes)
        v_pos_y = np.ones(v_nodes) * 0.8
        
        # Position check nodes at the bottom
        c_pos_x = np.linspace(0, 1, c_nodes)
        c_pos_y = np.ones(c_nodes) * 0.2
        
        # Plot variable nodes
        ax.scatter(v_pos_x, v_pos_y, s=100, c='blue', label='Variable Nodes')
        for i, (x, y) in enumerate(zip(v_pos_x, v_pos_y)):
            ax.text(x, y+0.05, f'v{i}', ha='center')
            
        # Plot check nodes
        ax.scatter(c_pos_x, c_pos_y, s=100, c='red', marker='s', label='Check Nodes')
        for i, (x, y) in enumerate(zip(c_pos_x, c_pos_y)):
            ax.text(x, y-0.05, f'c{i}', ha='center')
            
        # Plot connections
        for v, c in connections:
            ax.plot([v_pos_x[v], c_pos_x[c]], [v_pos_y[v], c_pos_y[c]], 'k-', alpha=0.3)
            
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(0, 1)
        ax.set_title('Tanner Graph Representation')
        ax.legend(loc='upper right')
        ax.axis('off')
        
    def plot_codeword(self, ax, codeword, title):
        """Plot a codeword/received word as a binary matrix."""
        ax.imshow(codeword, cmap='binary', aspect='auto')
        ax.set_title(title)
        ax.set_xlabel('Bit Position')
        ax.set_ylabel('Symbol Index')
        for i in range(codeword.shape[0]):
            for j in range(codeword.shape[1]):
                ax.text(j, i, str(codeword[i, j]), ha='center', va='center', color='r')
                
    def plot_llr_evolution(self, ax, llr_data, iteration):
        """Plot the evolution of LLR values for a specific iteration."""
        if iteration < len(llr_data):
            # Plot heatmap of LLR values
            im = ax.imshow(llr_data[iteration][:10, :10], cmap='coolwarm', aspect='auto')
            plt.colorbar(im, ax=ax)
            ax.set_title(f'LLR Values (Iteration {iteration+1})')
            ax.set_xlabel('GF Symbol')
            ax.set_ylabel('Code Symbol Index')
            
    def plot_convergence(self, ax):
        """Plot the convergence of the decoding process."""
        if self.iterations_data and self.fer_data and self.ber_data:
            # Real data from simulation
            ax.plot(self.frame_count_data, self.fer_data, 'b-', label='FER')
            ax.plot(self.frame_count_data, self.ber_data, 'r-', label='BER')
            ax.set_xlabel('Frame Count')
            ax.set_ylabel('Error Rate')
            ax.set_yscale('log')
            ax.legend()
            ax.set_title(f'Error Rates (SNR={self.snr_value}dB)')
            ax.grid(True, which='both', linestyle='--', alpha=0.5)
        else:
            # Synthetic data
            frames = np.arange(1, 101)
            fer = 0.5 * np.exp(-frames / 20)
            ber = 0.1 * np.exp(-frames / 30)
            
            ax.plot(frames, fer, 'b-', label='FER')
            ax.plot(frames, ber, 'r-', label='BER')
            ax.set_xlabel('Frame Count')
            ax.set_ylabel('Error Rate')
            ax.set_yscale('log')
            ax.legend()
            ax.set_title('Error Rates (Synthetic Data)')
            ax.grid(True, which='both', linestyle='--', alpha=0.5)
            
    def plot_iterations_histogram(self, ax):
        """Plot histogram of decoding iterations."""
        if len(self.iterations_data) > 0:
            # Real data
            avg_iter = np.mean(self.iterations_data)
            ax.hist(self.iterations_data, bins=10, alpha=0.7)
            ax.axvline(avg_iter, color='r', linestyle='--', label=f'Avg: {avg_iter:.2f}')
            ax.set_xlabel('Average Iterations')
            ax.set_ylabel('Frequency')
            ax.set_title('Decoding Iterations Distribution')
            ax.legend()
        else:
            # Synthetic data
            iterations = np.random.gamma(2, 1.5, 100)
            avg_iter = np.mean(iterations)
            ax.hist(iterations, bins=10, alpha=0.7)
            ax.axvline(avg_iter, color='r', linestyle='--', label=f'Avg: {avg_iter:.2f}')
            ax.set_xlabel('Average Iterations')
            ax.set_ylabel('Frequency')
            ax.set_title('Decoding Iterations Distribution')
            ax.legend()
            
    def plot_message_passing(self, ax, data, iteration):
        """Visualize message passing in the decoding process for a specific iteration."""
        # Create a simplified representation of message passing
        num_v = 8  # Number of variable nodes to show
        num_c = 4  # Number of check nodes to show
        
        # Position nodes
        v_pos_x = np.linspace(0.1, 0.9, num_v)
        v_pos_y = np.ones(num_v) * 0.8
        
        c_pos_x = np.linspace(0.2, 0.8, num_c)
        c_pos_y = np.ones(num_c) * 0.2
        
        # Plot nodes
        ax.scatter(v_pos_x, v_pos_y, s=100, c='blue', label='Variable Nodes')
        ax.scatter(c_pos_x, c_pos_y, s=100, c='red', marker='s', label='Check Nodes')
        
        # Show message passing with arrows
        # Direction depends on iteration phase (odd: V→C, even: C→V)
        if iteration % 2 == 0:
            # V→C messages
            for v in range(num_v):
                # Each variable node connects to 2 random check nodes
                for c in np.random.choice(num_c, size=2, replace=False):
                    ax.arrow(v_pos_x[v], v_pos_y[v]-0.05, 
                            c_pos_x[c]-v_pos_x[v], c_pos_y[c]-v_pos_y[v]+0.05, 
                            head_width=0.02, head_length=0.03, fc='green', ec='green', alpha=0.6)
            direction = "Variable to Check"
        else:
            # C→V messages
            for c in range(num_c):
                # Each check node connects to 2 random variable nodes
                for v in np.random.choice(num_v, size=2, replace=False):
                    ax.arrow(c_pos_x[c], c_pos_y[c]+0.05, 
                            v_pos_x[v]-c_pos_x[c], v_pos_y[v]-c_pos_y[c]-0.05, 
                            head_width=0.02, head_length=0.03, fc='purple', ec='purple', alpha=0.6)
            direction = "Check to Variable"
            
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(f'Message Passing - Iteration {iteration+1} ({direction})')
        ax.legend(loc='upper right')
        ax.axis('off')
        
    def plot_snr_vs_error(self, ax):
        """Plot error rates vs SNR (synthetic or from multiple runs)."""
        # Create synthetic data for SNR vs Error curve
        snr_values = np.linspace(0, 5, 6)
        fer_values = 0.5 * np.exp(-snr_values)
        ber_values = 0.1 * np.exp(-snr_values)
        
        ax.semilogy(snr_values, fer_values, 'bo-', label='FER')
        ax.semilogy(snr_values, ber_values, 'rx-', label='BER')
        
        # Mark current SNR point if available
        if self.snr_value > 0:
            if self.fer_data and self.ber_data:
                ax.plot(self.snr_value, self.fer_data[-1], 'bo', markersize=10, fillstyle='none')
                ax.plot(self.snr_value, self.ber_data[-1], 'rx', markersize=10, fillstyle='none')
                
        ax.set_xlabel('SNR (dB)')
        ax.set_ylabel('Error Rate')
        ax.set_title('Error Rate vs SNR')
        ax.grid(True, which='both', linestyle='--', alpha=0.5)
        ax.legend()
        
    def visualize_decoding_process(self, data=None):
        """Create a comprehensive visualization of the LDPC encoding, noise, and decoding process."""
        if data is None:
            data = self.generate_synthetic_data()
            
        iterations = data['iterations']
        
        # Create a figure with multiple subplots for different aspects
        fig = plt.figure(figsize=(15, 12))
        grid = plt.GridSpec(3, 3, figure=fig)
        
        # 1. Tanner graph representation
        ax_tanner = fig.add_subplot(grid[0, 0])
        self.plot_tanner_graph(ax_tanner, data['v_nodes'], data['c_nodes'], data['connections'])
        
        # 2. Original codeword
        ax_code = fig.add_subplot(grid[0, 1])
        self.plot_codeword(ax_code, data['codeword'], 'Original Codeword')
        
        # 3. Received (noisy) codeword
        ax_received = fig.add_subplot(grid[0, 2])
        self.plot_codeword(ax_received, data['received'], 'Received Word (with Noise)')
        
        # 4. LLR evolution for a specific iteration
        ax_llr = fig.add_subplot(grid[1, 0])
        self.plot_llr_evolution(ax_llr, data['llr_data'], 0)  # Show first iteration
        
        # 5. Message passing visualization
        ax_message = fig.add_subplot(grid[1, 1])
        self.plot_message_passing(ax_message, data, 0)  # Show first iteration
        
        # 6. Convergence graph
        ax_convergence = fig.add_subplot(grid[1, 2])
        self.plot_convergence(ax_convergence)
        
        # 7. Iterations histogram
        ax_iterations = fig.add_subplot(grid[2, 0])
        self.plot_iterations_histogram(ax_iterations)
        
        # 8. SNR vs Error curve
        ax_snr = fig.add_subplot(grid[2, 1:])
        self.plot_snr_vs_error(ax_snr)
        
        plt.tight_layout()
        
        # Add controls for iteration visualization
        axcolor = 'lightgoldenrodyellow'
        from matplotlib.widgets import Slider
        
        ax_slider = plt.axes([0.25, 0.02, 0.65, 0.03], facecolor=axcolor)
        slider = Slider(ax_slider, 'Iteration', 0, iterations-1, valinit=0, valstep=1)
        
        def update(val):
            iter_idx = int(slider.val)
            self.plot_llr_evolution(ax_llr, data['llr_data'], iter_idx)
            self.plot_message_passing(ax_message, data, iter_idx)
            fig.canvas.draw_idle()
            
        slider.on_changed(update)
        
        plt.suptitle('LDPC Code Simulation Visualization', fontsize=16)
        fig.tight_layout(rect=[0, 0.05, 1, 0.97])
        
        return fig
        
    def run_interactive_dashboard(self):
        """Run an interactive dashboard for the LDPC simulation."""
        import ipywidgets as widgets
        from IPython.display import display
        
        # Generate synthetic data
        data = self.generate_synthetic_data()
        
        # Create widgets
        iter_slider = widgets.IntSlider(
            value=0,
            min=0,
            max=data['iterations']-1,
            step=1,
            description='Iteration:',
            continuous_update=False
        )
        
        output = widgets.Output()
        
        # Create tabs for different visualizations
        tab1 = widgets.Output()
        tab2 = widgets.Output()
        tab3 = widgets.Output()
        
        tabs = widgets.Tab(children=[tab1, tab2, tab3])
        tabs.set_title(0, 'Code Structure')
        tabs.set_title(1, 'Decoding Process')
        tabs.set_title(2, 'Performance')
        
        def on_iter_change(change):
            with output:
                output.clear_output(wait=True)
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
                self.plot_llr_evolution(ax1, data['llr_data'], change['new'])
                self.plot_message_passing(ax2, data, change['new'])
                plt.tight_layout()
                plt.show()
                
        iter_slider.observe(on_iter_change, names='value')
        
        # Create content for tab 1 (Code Structure)
        with tab1:
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
            self.plot_tanner_graph(ax1, data['v_nodes'], data['c_nodes'], data['connections'])
            self.plot_codeword(ax2, data['codeword'], 'Original Codeword')
            self.plot_codeword(ax3, data['received'], 'Received Word (with Noise)')
            plt.tight_layout()
            plt.show()
            
        # Create content for tab 2 (Decoding Process)
        with tab2:
            display(iter_slider)
            display(output)
            on_iter_change({'new': 0})  # Initialize with iteration 0
            
        # Create content for tab 3 (Performance)
        with tab3:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            self.plot_convergence(ax1)
            self.plot_snr_vs_error(ax2)
            plt.tight_layout()
            plt.show()
            
        display(tabs)

    def create_animation(self, output_file='ldpc_decoding.mp4'):
        """Create an animation of the decoding process."""
        import matplotlib.animation as animation
        
        data = self.generate_synthetic_data()
        iterations = data['iterations']
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        def update_frame(i):
            ax1.clear()
            ax2.clear()
            
            self.plot_llr_evolution(ax1, data['llr_data'], i % iterations)
            self.plot_message_passing(ax2, data, i % iterations)
            
            return ax1, ax2
            
        ani = animation.FuncAnimation(fig, update_frame, frames=iterations, interval=1000, blit=False)
        
        # Save animation
        ani.save(output_file, writer='ffmpeg', fps=1)
        
        return ani

    def extract_real_data_from_c_output(self, output_file):
        """Extract real data from C program output for direct visualization."""
        try:
            with open(output_file, 'r') as file:
                content = file.read()
                self.parse_c_output(content)
                return True
        except Exception as e:
            print(f"Error reading output file: {e}")
            return False

def main():
    """Main function to execute the LDPC visualization."""
    visualizer = LDPCVisualizer()
    
    # You can extract data from a saved output file if you have one
    # visualizer.extract_real_data_from_c_output("simulation_output.txt")
    
    # Or run the C program directly (uncomment and adjust path as needed)
    # visualizer.run_simulation("./ldpc_decoder", ["100", "10", "matrix.txt", "2.0", "3", "3", "0.5", "100"])
    
    # For demonstration, use synthetic data
    fig = visualizer.visualize_decoding_process()
    
    # Save the figure
    fig.savefig("./data/ldpc_visualization.png", dpi=300, bbox_inches='tight')
    
    # For interactive use in Jupyter notebook, uncomment:
    # visualizer.run_interactive_dashboard()
    
    # Create animation
    # visualizer.create_animation()
    
    # Show the plot
    plt.show()
    
    print("Visualization complete!")

if __name__ == "__main__":
    main()