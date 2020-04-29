Imports System.Drawing
Imports System.Drawing.Drawing2D
Imports System.Drawing.Imaging
Imports System.Drawing.Text

'Public Class Win32
'    Public Declare Function BitBlt Lib "gdi32" Alias "BitBlt" (ByVal hDestDC As Integer, ByVal x As Integer, ByVal y As Integer, ByVal nWidth As Integer, ByVal nHeight As Integer, ByVal hSrcDC As Integer, ByVal xSrc As Integer, ByVal ySrc As Integer, ByVal dwRop As Integer) As Integer

'    Public Declare Function GetWindowDC Lib "user32" Alias "GetWindowDC" (ByVal hwnd As Integer) As Integer

'    Public Declare Function ReleaseDC Lib "user32" Alias "ReleaseDC" (ByVal hwnd As Integer, ByVal hdc As Integer) As Integer

'    Public Const SRCCOPY As Integer = &HCC0020
'End Class

Public Class frmChlor

    Inherits System.Windows.Forms.Form

#Region " Windows Form Designer generated code "

    Public Sub New()
        MyBase.New()

        'This call is required by the Windows Form Designer.
        InitializeComponent()

        'Add any initialization after the InitializeComponent() call

    End Sub

    'Form overrides dispose to clean up the component list.
    Protected Overloads Overrides Sub Dispose(ByVal disposing As Boolean)
        If disposing Then
            If Not (components Is Nothing) Then
                components.Dispose()
            End If
        End If
        MyBase.Dispose(disposing)
    End Sub

    'Required by the Windows Form Designer
    Private components As System.ComponentModel.IContainer

    'NOTE: The following procedure is required by the Windows Form Designer
    'It can be modified using the Windows Form Designer.  
    'Do not modify it using the code editor.
    Public WithEvents Command1 As System.Windows.Forms.Button
    Public WithEvents ToolTip1 As System.Windows.Forms.ToolTip
    Friend WithEvents PictureBox1 As System.Windows.Forms.PictureBox
    Friend WithEvents PictureBox2 As System.Windows.Forms.PictureBox
    Friend WithEvents PictureBox3 As System.Windows.Forms.PictureBox
    Friend WithEvents PictureBox4 As System.Windows.Forms.PictureBox
    Friend WithEvents PictureBox5 As System.Windows.Forms.PictureBox
    Friend WithEvents PictureBox6 As System.Windows.Forms.PictureBox
    <System.Diagnostics.DebuggerStepThrough()> Private Sub InitializeComponent()
        Me.components = New System.ComponentModel.Container
        Me.Command1 = New System.Windows.Forms.Button
        Me.ToolTip1 = New System.Windows.Forms.ToolTip(Me.components)
        Me.PictureBox1 = New System.Windows.Forms.PictureBox
        Me.PictureBox2 = New System.Windows.Forms.PictureBox
        Me.PictureBox3 = New System.Windows.Forms.PictureBox
        Me.PictureBox4 = New System.Windows.Forms.PictureBox
        Me.PictureBox5 = New System.Windows.Forms.PictureBox
        Me.PictureBox6 = New System.Windows.Forms.PictureBox
        Me.SuspendLayout()
        '
        'Command1
        '
        Me.Command1.BackColor = System.Drawing.SystemColors.Control
        Me.Command1.Cursor = System.Windows.Forms.Cursors.Default
        Me.Command1.Font = New System.Drawing.Font("Arial", 8.0!, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, CType(0, Byte))
        Me.Command1.ForeColor = System.Drawing.SystemColors.ControlText
        Me.Command1.Location = New System.Drawing.Point(872, 456)
        Me.Command1.Name = "Command1"
        Me.Command1.RightToLeft = System.Windows.Forms.RightToLeft.No
        Me.Command1.Size = New System.Drawing.Size(89, 25)
        Me.Command1.TabIndex = 19
        Me.Command1.Text = "Store"
        '
        'PictureBox1
        '
        Me.PictureBox1.BackColor = System.Drawing.Color.White
        Me.PictureBox1.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle
        Me.PictureBox1.Location = New System.Drawing.Point(248, 88)
        Me.PictureBox1.Name = "PictureBox1"
        Me.PictureBox1.TabIndex = 20
        Me.PictureBox1.TabStop = False
        '
        'PictureBox2
        '
        Me.PictureBox2.BackColor = System.Drawing.Color.White
        Me.PictureBox2.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle
        Me.PictureBox2.Location = New System.Drawing.Point(360, 88)
        Me.PictureBox2.Name = "PictureBox2"
        Me.PictureBox2.TabIndex = 21
        Me.PictureBox2.TabStop = False
        '
        'PictureBox3
        '
        Me.PictureBox3.BackColor = System.Drawing.Color.White
        Me.PictureBox3.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle
        Me.PictureBox3.Location = New System.Drawing.Point(472, 88)
        Me.PictureBox3.Name = "PictureBox3"
        Me.PictureBox3.TabIndex = 22
        Me.PictureBox3.TabStop = False
        '
        'PictureBox4
        '
        Me.PictureBox4.BackColor = System.Drawing.Color.White
        Me.PictureBox4.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle
        Me.PictureBox4.Location = New System.Drawing.Point(584, 88)
        Me.PictureBox4.Name = "PictureBox4"
        Me.PictureBox4.TabIndex = 23
        Me.PictureBox4.TabStop = False
        '
        'PictureBox5
        '
        Me.PictureBox5.BackColor = System.Drawing.Color.White
        Me.PictureBox5.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle
        Me.PictureBox5.Location = New System.Drawing.Point(696, 88)
        Me.PictureBox5.Name = "PictureBox5"
        Me.PictureBox5.TabIndex = 24
        Me.PictureBox5.TabStop = False
        '
        'PictureBox6
        '
        Me.PictureBox6.BackColor = System.Drawing.Color.White
        Me.PictureBox6.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle
        Me.PictureBox6.Location = New System.Drawing.Point(808, 88)
        Me.PictureBox6.Name = "PictureBox6"
        Me.PictureBox6.TabIndex = 25
        Me.PictureBox6.TabStop = False
        '
        'frmChlor
        '
        Me.AutoScaleBaseSize = New System.Drawing.Size(5, 13)
        Me.ClientSize = New System.Drawing.Size(968, 485)
        Me.Controls.Add(Me.PictureBox1)
        Me.Controls.Add(Me.PictureBox5)
        Me.Controls.Add(Me.PictureBox4)
        Me.Controls.Add(Me.PictureBox3)
        Me.Controls.Add(Me.PictureBox2)
        Me.Controls.Add(Me.Command1)
        Me.Controls.Add(Me.PictureBox6)
        Me.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle
        Me.Name = "frmChlor"
        Me.Text = "frmChlor"
        Me.WindowState = System.Windows.Forms.FormWindowState.Maximized
        Me.ResumeLayout(False)

    End Sub

#End Region

    Dim a As Single
    Dim Picture_Number As Short = 0

    'lors du chargement de la fenêtre
    Private Sub frmChlor_Load(ByVal sender As Object, ByVal e As System.EventArgs) Handles MyBase.Load
        a = Me.Width / 34
        Dim x As Short = a
        Dim y As Short = 18 * a
        Me.Command1.Location = New System.Drawing.Point(x, y)
    End Sub

    'dessin des graphiques durant l'exécution du calcul
    Public Sub design(ByRef wsat As Single, ByRef H_old() As Decimal, ByRef W() As Decimal, ByRef T_old() As Decimal, ByRef C_new() As Decimal, ByRef CT() As Decimal, ByRef tempmin As Single, ByRef tempmax As Single, ByRef Length As Single, ByRef PosProf() As Decimal, ByRef hMin As Single, ByRef hEcart As Single, ByRef wMin As Single, ByRef wEcart As Single, ByRef CTmin As Single, ByRef CTecart As Single, ByRef CTmax As Integer, ByRef CLmin As Single, ByRef CLecart As Single, ByRef CLmax As Integer, ByRef Tecart As Single, ByRef Dofs As Short, ByRef Tijd As Decimal, ByRef Gxc As Single, ByRef Dxc As Single, ByRef Ph() As Decimal)
        ' Get a Graphics object from the current form and clear its background.
        Dim x1 As Single
        Dim x2 As Single
        Dim h1 As Single
        Dim h2 As Single
        Dim i As Short
        Dim EctX As Single
        Dim EctY As Single
        Dim Tmax As Single
        Dim m_Longueur As Single
        Dim m_Hauteur As Single

        Dim position As Point
        Dim couleur As Color
        Dim Message0 As String
        Dim Message1 As String
        Dim Message2 As String
        Dim Message3 As String

        a = Me.Width / 34

        ''''''''''''''''''''''''''''''''''''
        ' Construction du premier graphique
        ''''''''''''''''''''''''''''''''''''
        position = New Point(a, a)
        m_Longueur = a * 10             'ne change pas pour les autres graphiques
        m_Hauteur = a * 7.5             'ne change pas pour les autres graphiques

        couleur = System.Drawing.Color.Red
        Message0 = ""                       'ne change pas pour les autres graphiques
        Message1 = "Temperature Potential [°C]"
        Message2 = "Depth [mm]"             'ne change pas pour les autres grapiques
        Message3 = "Temperature Potential Distribution at " & Format(Tijd, "0.00") & " Days"
        x1 = 0                              'ne change pas pour les autres graphiques
        x2 = Length                         'ne change pas pour les autres graphiques
        h1 = (Int(tempmin / 10)) * 10
        If Int(tempmax / 10) - tempmax / 10 = 0 Then
            Tmax = tempmax
        Else
            Tmax = (Int(tempmax / 10) + 1) * 10
        End If
        h2 = Tmax
        EctX = Int(Length / 10)             'ne change pas pour les autres graphiques
        EctY = Tecart
        Dim gr01 As CGraph = New CGraph(Me, PictureBox1, position, m_Longueur, m_Hauteur, couleur, Message0, Message1, Message2, Message3, x1, x2, h1, h2, EctX, EctY, PosProf, T_old, Dofs, -1, -1)
        ''''''''''''''''''''''''''''''''''''
        ' Construction du deuxième graphique
        ''''''''''''''''''''''''''''''''''''
        position = New Point(12 * a, a)
        couleur = System.Drawing.Color.BlueViolet
        Message1 = "Moisture Potential [P/Ps]"
        Message3 = "Moisture Potential Distribution at " & Format(Tijd, "0.00") & " Days"
        h1 = hMin
        h2 = 1
        EctY = 0.1
        Dim gr02 As CGraph = New CGraph(Me, PictureBox2, position, m_Longueur, m_Hauteur, couleur, Message0, Message1, Message2, Message3, x1, x2, h1, h2, EctX, EctY, PosProf, H_old, Dofs, -1, -1)
        ''''''''''''''''''''''''''''''''''''
        ' Construction du troisième graphique
        ''''''''''''''''''''''''''''''''''''
        position = New Point(12 * a, 10 * a)
        couleur = System.Drawing.Color.Blue
        Message1 = "Moisture Content [kg/m3]"
        Message3 = "Moisture Content Distribution at " & Format(Tijd, "0.00") & " Days"
        h1 = wMin
        i = CShort(wsat / 10.0#)
        If wsat > CSng(i) * CSng(10) Then
            i = i + CShort(1)
        End If
        h2 = i * 10
        EctY = wEcart
        Dim gr03 As CGraph = New CGraph(Me, PictureBox3, position, m_Longueur, m_Hauteur, couleur, Message0, Message1, Message2, Message3, x1, x2, h1, h2, EctX, EctY, PosProf, W, Dofs, -1, -1)
        ''''''''''''''''''''''''''''''''''''
        ' Construction du quatrième graphique
        ''''''''''''''''''''''''''''''''''''
        position = New Point(23 * a, a)
        couleur = System.Drawing.Color.GreenYellow
        Message1 = "Total Chloride Ion Content [kg/m3]"
        Message3 = "Total Chloride Ion Distribution at " & Format(Tijd, "0.00") & " Days"
        h1 = CTmin
        h2 = CTmax
        If (h2 - h1) / CTecart > 14 Or (h2 - h1) / CTecart < 2 Then
            EctY = (h2 - h1) / 10
        Else
            EctY = CTecart
        End If
        Dim gr04 As CGraph = New CGraph(Me, PictureBox4, position, m_Longueur, m_Hauteur, couleur, Message0, Message1, Message2, Message3, x1, x2, h1, h2, EctX, EctY, PosProf, CT, Dofs, -1, -1)
        ''''''''''''''''''''''''''''''''''''
        ' Construction du cinquième graphique
        ''''''''''''''''''''''''''''''''''''
        position = New Point(23 * a, 10 * a)
        couleur = System.Drawing.Color.Green
        Message1 = "Free Chloride Ion Content [kg/m3]"
        Message3 = "Free Chloride Ion Distribution at " & Format(Tijd, "0.00") & " Days"
        h1 = CLmin
        h2 = CLmax
        If (h2 - h1) / CLecart > 14 Or (h2 - h1) / CLecart < 2 Then
            EctY = (h2 - h1) / 10
        Else
            EctY = CLecart
        End If
        Dim DTot(Dofs) As Decimal
        For i = 1 To Dofs
            DTot(i) = C_new(i) * W(i) / 1000.0#
        Next i
        Dim gr05 As CGraph = New CGraph(Me, PictureBox5, position, m_Longueur, m_Hauteur, couleur, Message0, Message1, Message2, Message3, x1, x2, h1, h2, EctX, EctY, PosProf, DTot, Dofs, -1, -1)
        ''''''''''''''''''''''''''''''''''''
        ' Constuction du sixième graphique
        ''''''''''''''''''''''''''''''''''''
        position = New Point(a, 10 * a)

        couleur = System.Drawing.Color.DarkCyan
        Message1 = "PH"
        Message3 = "PH Distribution at " & Format(Tijd, "0.00") & " Days"
        h1 = 0
        h2 = 14
        EctY = 1
        Dim gr06 As CGraph = New CGraph(Me, PictureBox6, position, m_Longueur, m_Hauteur, couleur, Message0, Message1, Message2, Message3, x1, x2, h1, h2, EctX, EctY, PosProf, Ph, Dofs, Gxc, Dxc)
    End Sub

    'Click sur le bouton store
    Private Sub Command1_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Command1.Click
        Save_pictures(Picture_Number)
    End Sub

End Class
