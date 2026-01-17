import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, Button
from mpl_toolkits.mplot3d import Axes3D

class ElectricFieldViz:
    def __init__(self):
        # Dark theme setup
        self.fig = plt.figure(figsize=(12, 8), facecolor='black')
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_facecolor('black')
        plt.subplots_adjust(left=0.25, bottom=0.25)
        
        # Initial Physics Parameters
        self.x0 = 2.0  # Observer starts on X axis
        self.R0 = 1.0
        self.Q0 = 1.0
        self.N0 = 20
        
        # State
        self.is_disk = False
        self.show_vectors = True
        self.show_net = True
        self.follow_observer = False
        
        # Camera State
        # Re-oriented default view: Side view showing the ring profile
        self.azim = -60.0 
        self.elev = 20.0
        self.dist = 4.0      
        self.target = np.array([0.0, 0.0, 0.0])
        self.mouse_drag_active = False
        self.last_mouse_x = 0
        self.last_mouse_y = 0
        self.mouse_btn = None 
        
        # UI Setup
        self.setup_widgets()
        
        # Event binding
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        
        # Initial Draw
        self.ax.view_init(elev=self.elev, azim=self.azim)
        self.update_geometry()
        
    def setup_widgets(self):
        # Style props
        slider_bg = '#404040'
        text_color = 'white'
        
        # Sliders
        # Changed Z slider to X slider
        ax_x = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=slider_bg)
        ax_R = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor=slider_bg)
        ax_Q = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=slider_bg)
        ax_N = plt.axes([0.05, 0.25, 0.02, 0.63], facecolor=slider_bg)
        
        self.slider_x = Slider(ax_x, 'X Position', -5.0, 5.0, valinit=self.x0, color='gray')
        self.slider_R = Slider(ax_R, 'Radius R', 0.5, 3.0, valinit=self.R0, color='gray')
        self.slider_Q = Slider(ax_Q, 'Charge Q', -5.0, 5.0, valinit=self.Q0, color='gray')
        self.slider_N = Slider(ax_N, 'Count N', 4, 100, valinit=self.N0, valstep=1, orientation='vertical', color='gray')
        
        # Styling texts
        for s in [self.slider_x, self.slider_R, self.slider_Q, self.slider_N]:
            s.label.set_color(text_color)
            s.valtext.set_color(text_color)
        
        self.slider_x.on_changed(self.on_param_change)
        self.slider_R.on_changed(self.on_param_change)
        self.slider_Q.on_changed(self.on_param_change)
        self.slider_N.on_changed(self.on_param_change)
        
        # Checkboxes
        # Top-Right. frameon=True with light color to ensures text is visible 
        # regardless of matplotlib version/styling hacks.
        ax_check = plt.axes([0.80, 0.80, 0.15, 0.15], frameon=True, facecolor='#dddddd')
        self.check = CheckButtons(ax_check, ['Mode: Disk', 'Show Vectors', 'Show Net', 'Follow Obs'], [False, True, True, False])
        
        # No more fragile styling hacks. Default black text works on light gray.
            
        self.check.on_clicked(self.update_check)
        
        # Reset View Button
        ax_reset = plt.axes([0.8, 0.02, 0.1, 0.04])
        self.btn_reset = Button(ax_reset, 'Reset View', color='#404040', hovercolor='#606060')
        self.btn_reset.label.set_color('white')
        self.btn_reset.on_clicked(self.reset_view)
        
    def reset_view(self, event):
        self.azim = -60.0
        self.elev = 20.0
        self.dist = 4.0
        self.target = np.array([0.0, 0.0, 0.0])
        self.update_camera()

    def update_check(self, label):
        if label == 'Mode: Disk':
            self.is_disk = not self.is_disk
        elif label == 'Show Vectors':
            self.show_vectors = not self.show_vectors
        elif label == 'Show Net':
            self.show_net = not self.show_net
        elif label == 'Follow Obs':
            self.follow_observer = not self.follow_observer
        self.update_geometry()

    def on_param_change(self, val):
        self.update_geometry()

    # --- Interaction Events ---
    def on_press(self, event):
        if event.inaxes != self.ax: return
        self.mouse_drag_active = True
        self.last_mouse_x = event.x
        self.last_mouse_y = event.y
        self.mouse_btn = event.button

    def on_release(self, event):
        self.mouse_drag_active = False
        self.mouse_btn = None

    def on_motion(self, event):
        if not self.mouse_drag_active or event.inaxes != self.ax:
            return
            
        dx = event.x - self.last_mouse_x
        dy = event.y - self.last_mouse_y
        self.last_mouse_x = event.x
        self.last_mouse_y = event.y
        
        if self.mouse_btn == 1: # Left: Orbit
            sensitivity = 0.5
            self.azim -= dx * sensitivity
            self.elev += dy * sensitivity 
            self.elev = np.clip(self.elev, -90, 90)
            self.update_camera()
            
        elif self.mouse_btn == 3: # Right: Pan
            pan_sens = self.dist * 0.005
            self.target[1] -= dx * pan_sens
            self.target[2] -= dy * pan_sens
            self.update_camera()

    def on_scroll(self, event):
        if event.inaxes != self.ax: return
        zoom_sensitivity = 0.1
        if event.step > 0:
            self.dist *= (1.0 - zoom_sensitivity)
        else:
            self.dist *= (1.0 + zoom_sensitivity)
        self.update_camera()

    def update_camera(self):
        self.ax.view_init(elev=self.elev, azim=self.azim)
        
        t = self.target
        d = self.dist
        
        # Enforce "infinite space" vibe by tight clipping around target
        self.ax.set_xlim(t[0] - d, t[0] + d)
        self.ax.set_ylim(t[1] - d, t[1] + d)
        self.ax.set_zlim(t[2] - d, t[2] + d)
        self.ax.set_box_aspect([1, 1, 1])
        self.fig.canvas.draw_idle()

    def generate_charge_distribution(self, N, R):
        if not self.is_disk:
            # Ring in Y-Z plane
            theta = np.linspace(0, 2*np.pi, int(N), endpoint=False)
            x = np.zeros_like(theta)
            y = R * np.cos(theta)
            z = R * np.sin(theta)
            return x, y, z
        else:
            # Disk in Y-Z plane
            indices = np.arange(0, int(N), dtype=float) + 0.5
            r = R * np.sqrt(indices/int(N))
            theta = np.pi * (1 + 5**0.5) * indices
            
            x = np.zeros_like(r)
            y = r * np.cos(theta)
            z = r * np.sin(theta)
            return x, y, z

    def update_geometry(self):
        x_obs = self.slider_x.val
        R = self.slider_R.val
        Q = self.slider_Q.val
        N = int(self.slider_N.val)
        
        # Update target if following
        if self.follow_observer:
            # Follow X
            self.target[0] = x_obs
            
        self.ax.clear()
        
        # Infinite space styling (re-apply after clear)
        self.ax.set_axis_off() 
        self.ax.grid(False)
        self.ax.set_facecolor('black')
        
        # Draw minimal custom axes at origin
        len_ax = 3.0
        self.ax.plot([-len_ax, len_ax], [0, 0], [0, 0], 'w-', lw=0.5, alpha=0.3) # X
        self.ax.plot([0, 0], [-len_ax, len_ax], [0, 0], 'w-', lw=0.5, alpha=0.3) # Y
        self.ax.plot([0, 0], [0, 0], [-len_ax, len_ax], 'w-', lw=0.5, alpha=0.3) # Z
        self.ax.text(len_ax, 0, 0, "x", color='white', fontsize=8)
        self.ax.text(0, len_ax, 0, "y", color='white', fontsize=8)
        self.ax.text(0, 0, len_ax, "z", color='white', fontsize=8)
        
        # Observer: Now on X axis
        self.ax.scatter([x_obs], [0], [0], color='#ff3333', s=60, label='Observer', zorder=100, edgecolors='white', linewidth=0.5)
        
        # Charges
        xc, yc, zc = self.generate_charge_distribution(N, R)
        if self.is_disk:
            self.ax.scatter(xc, yc, zc, c='cyan', alpha=0.6, s=15, edgecolors='none')
            # Rim in Y-Z
            t = np.linspace(0, 2*np.pi, 60)
            self.ax.plot(np.zeros_like(t), R*np.cos(t), R*np.sin(t), c='cyan', alpha=0.4, lw=1)
        else:
            self.ax.scatter(xc, yc, zc, c='cyan', s=25, edgecolors='white', linewidth=0.5)
            # Wireframe ring in Y-Z
            t = np.linspace(0, 2*np.pi, 60)
            self.ax.plot(np.zeros_like(t), R*np.cos(t), R*np.sin(t), c='cyan', lw=2)

        # Fields
        dQ = Q / N
        P = np.array([x_obs, 0, 0])
        Ex_net, Ey_net, Ez_net = 0, 0, 0
        
        # Separate Scales
        scale_dE = 50.0 # Huge for components
        scale_net = 10.0 # Normal for net
        
        if self.show_vectors:
            for i in range(N):
                source = np.array([xc[i], yc[i], zc[i]])
                r_vec = P - source
                r_mag = np.linalg.norm(r_vec)
                if r_mag < 1e-6: continue 
                
                r_hat = r_vec / r_mag
                dE_mag = abs(dQ) / (r_mag**2)
                
                if Q < 0:
                    dE_vec = -dE_mag * r_hat
                else:
                    dE_vec = dE_mag * r_hat
                
                Ex_net += dE_vec[0]
                Ey_net += dE_vec[1]
                Ez_net += dE_vec[2]
                
                # Visual vector
                # High contrast Yellow, thicker, opaque
                self.ax.quiver(P[0], P[1], P[2], 
                               dE_vec[0], dE_vec[1], dE_vec[2],
                               length=np.linalg.norm(dE_vec)*scale_dE, normalize=True,
                               color='yellow', alpha=0.5, arrow_length_ratio=0.1, linewidth=1.5)
        else:
            for i in range(N):
                source = np.array([xc[i], yc[i], zc[i]])
                r_vec = P - source
                r_mag = np.linalg.norm(r_vec)
                if r_mag < 1e-6: continue
                r_hat = r_vec / r_mag
                if Q < 0:
                    dE_vec = -(abs(dQ) / (r_mag**2)) * r_hat
                else:
                    dE_vec = (abs(dQ) / (r_mag**2)) * r_hat
                Ex_net += dE_vec[0]
                Ey_net += dE_vec[1]
                Ez_net += dE_vec[2]

        if self.show_net:
            net_mag = np.sqrt(Ex_net**2 + Ey_net**2 + Ez_net**2)
            if net_mag > 1e-9:
                self.ax.quiver(P[0], P[1], P[2],
                               Ex_net, Ey_net, Ez_net,
                               length=net_mag*scale_net, normalize=True,
                               color='#00ff00', linewidth=4.0, arrow_length_ratio=0.2)
                               
        self.update_camera()

if __name__ == "__main__":
    app = ElectricFieldViz()
    plt.show()
