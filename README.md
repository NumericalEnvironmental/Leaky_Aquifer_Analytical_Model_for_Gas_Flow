# Leaky_Aquifer_Analytical_Model_for_Gas_Flow

![Preview](https://numericalenvironmental.files.wordpress.com/2020/04/flow-field.jpg?w=1632)

This Python script models the steady-state pressure distribution and gas mass fluxes associated with operation of a soil vapor extraction (SVE) system in a leaky aquifer setting. The analytical solution is based on a radial groundwater flow model that has been modified for application to a compressible gas; the details of the model and an example application are provided in a blog post https://github.com/NumericalEnvironmental/Leaky_Aquifer_Analytical_Model_for_Gas_Flow. The script also uses SciPy’s multivariable nonlinear optimizer (https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.fsolve.html) and numerical integration routine (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html) to simulate passive air intake wells and to computer gas fluxes normal to an encompassing polygon, respectively.
The following text-based input files are required by the script:
* reservoir.txt – aquifer and aquitard thicknesses and permeabilities
* gas.txt – gas-phase properties (e.g., viscosity, molecular weight)
* wells.csv – locations and extraction rates for SVE wells
* air intake wells.txt – locations and (placeholder) injection rates for passive air intake wells
* polygon.csv – points defining a bounding polygon encompassing the SVE and air intake well arrays
* grid.txt – gridding parameters for setting up model output (for subsequent contouring/plotting)

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
