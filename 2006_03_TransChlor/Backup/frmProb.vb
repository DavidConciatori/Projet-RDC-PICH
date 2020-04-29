Public Class frmProb
    Inherits System.Windows.Forms.Form
    Dim Length As Single
    Dim Labe(1) As Object
    Dim Lambda(1, 1, 1) As Double
    Dim Ksi(1, 1, 1) As Double
    Dim Nline As Short
    Dim Nc As Short
    Dim PreFile(7) As String
    Dim OutFile As String
    Dim Files As Short = 0

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
    Friend WithEvents GroupBox1 As System.Windows.Forms.GroupBox
    Friend WithEvents GroupBox3 As System.Windows.Forms.GroupBox
    Friend WithEvents GroupBox2 As System.Windows.Forms.GroupBox
    Friend WithEvents CheckBox6 As System.Windows.Forms.CheckBox
    Friend WithEvents CheckBox7 As System.Windows.Forms.CheckBox
    Friend WithEvents CheckBox8 As System.Windows.Forms.CheckBox
    Friend WithEvents Label1 As System.Windows.Forms.Label
    Friend WithEvents Label2 As System.Windows.Forms.Label
    Friend WithEvents Label3 As System.Windows.Forms.Label
    Friend WithEvents Label4 As System.Windows.Forms.Label
    Friend WithEvents TextBox1 As System.Windows.Forms.TextBox
    Friend WithEvents TextBox2 As System.Windows.Forms.TextBox
    Friend WithEvents ComboBox1 As System.Windows.Forms.ComboBox
    Friend WithEvents Label5 As System.Windows.Forms.Label
    Friend WithEvents Label6 As System.Windows.Forms.Label
    Friend WithEvents Label7 As System.Windows.Forms.Label
    Friend WithEvents Label8 As System.Windows.Forms.Label
    Friend WithEvents Label9 As System.Windows.Forms.Label
    Friend WithEvents Button1 As System.Windows.Forms.Button
    Friend WithEvents Button2 As System.Windows.Forms.Button
    Friend WithEvents TextBox3 As System.Windows.Forms.TextBox
    Friend WithEvents TextBox4 As System.Windows.Forms.TextBox
    Friend WithEvents Button3 As System.Windows.Forms.Button
    Friend WithEvents Label10 As System.Windows.Forms.Label
    Friend WithEvents Label11 As System.Windows.Forms.Label
    Friend WithEvents Label12 As System.Windows.Forms.Label
    Friend WithEvents Label13 As System.Windows.Forms.Label
    Friend WithEvents CheckBox5 As System.Windows.Forms.CheckBox
    Friend WithEvents CheckBox4 As System.Windows.Forms.CheckBox
    Friend WithEvents CheckBox3 As System.Windows.Forms.CheckBox
    Friend WithEvents CheckBox2 As System.Windows.Forms.CheckBox
    Friend WithEvents CheckBox1 As System.Windows.Forms.CheckBox
    <System.Diagnostics.DebuggerStepThrough()> Private Sub InitializeComponent()
        Me.GroupBox1 = New System.Windows.Forms.GroupBox
        Me.CheckBox7 = New System.Windows.Forms.CheckBox
        Me.CheckBox6 = New System.Windows.Forms.CheckBox
        Me.CheckBox5 = New System.Windows.Forms.CheckBox
        Me.CheckBox4 = New System.Windows.Forms.CheckBox
        Me.CheckBox3 = New System.Windows.Forms.CheckBox
        Me.CheckBox2 = New System.Windows.Forms.CheckBox
        Me.CheckBox1 = New System.Windows.Forms.CheckBox
        Me.GroupBox3 = New System.Windows.Forms.GroupBox
        Me.Label13 = New System.Windows.Forms.Label
        Me.Label12 = New System.Windows.Forms.Label
        Me.Label11 = New System.Windows.Forms.Label
        Me.Label10 = New System.Windows.Forms.Label
        Me.TextBox4 = New System.Windows.Forms.TextBox
        Me.TextBox3 = New System.Windows.Forms.TextBox
        Me.Label9 = New System.Windows.Forms.Label
        Me.Label8 = New System.Windows.Forms.Label
        Me.Label7 = New System.Windows.Forms.Label
        Me.Label6 = New System.Windows.Forms.Label
        Me.Label5 = New System.Windows.Forms.Label
        Me.ComboBox1 = New System.Windows.Forms.ComboBox
        Me.GroupBox2 = New System.Windows.Forms.GroupBox
        Me.TextBox2 = New System.Windows.Forms.TextBox
        Me.TextBox1 = New System.Windows.Forms.TextBox
        Me.Label4 = New System.Windows.Forms.Label
        Me.Label3 = New System.Windows.Forms.Label
        Me.Label2 = New System.Windows.Forms.Label
        Me.Label1 = New System.Windows.Forms.Label
        Me.CheckBox8 = New System.Windows.Forms.CheckBox
        Me.Button1 = New System.Windows.Forms.Button
        Me.Button2 = New System.Windows.Forms.Button
        Me.Button3 = New System.Windows.Forms.Button
        Me.GroupBox1.SuspendLayout()
        Me.GroupBox3.SuspendLayout()
        Me.GroupBox2.SuspendLayout()
        Me.SuspendLayout()
        '
        'GroupBox1
        '
        Me.GroupBox1.Controls.Add(Me.CheckBox7)
        Me.GroupBox1.Controls.Add(Me.CheckBox6)
        Me.GroupBox1.Controls.Add(Me.CheckBox5)
        Me.GroupBox1.Controls.Add(Me.CheckBox4)
        Me.GroupBox1.Controls.Add(Me.CheckBox3)
        Me.GroupBox1.Controls.Add(Me.CheckBox2)
        Me.GroupBox1.Controls.Add(Me.CheckBox1)
        Me.GroupBox1.Location = New System.Drawing.Point(8, 8)
        Me.GroupBox1.Name = "GroupBox1"
        Me.GroupBox1.Size = New System.Drawing.Size(184, 248)
        Me.GroupBox1.TabIndex = 0
        Me.GroupBox1.TabStop = False
        Me.GroupBox1.Text = "Traitement probabiliste"
        '
        'CheckBox7
        '
        Me.CheckBox7.Checked = True
        Me.CheckBox7.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox7.Location = New System.Drawing.Point(8, 216)
        Me.CheckBox7.Name = "CheckBox7"
        Me.CheckBox7.Size = New System.Drawing.Size(112, 24)
        Me.CheckBox7.TabIndex = 6
        Me.CheckBox7.Text = "Chlorures totaux"
        '
        'CheckBox6
        '
        Me.CheckBox6.Checked = True
        Me.CheckBox6.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox6.Location = New System.Drawing.Point(8, 184)
        Me.CheckBox6.Name = "CheckBox6"
        Me.CheckBox6.TabIndex = 5
        Me.CheckBox6.Text = "Chlorures libres"
        '
        'CheckBox5
        '
        Me.CheckBox5.Checked = True
        Me.CheckBox5.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox5.Location = New System.Drawing.Point(8, 152)
        Me.CheckBox5.Name = "CheckBox5"
        Me.CheckBox5.Size = New System.Drawing.Size(168, 24)
        Me.CheckBox5.TabIndex = 4
        Me.CheckBox5.Text = "Profondeur de carbonatation"
        '
        'CheckBox4
        '
        Me.CheckBox4.Checked = True
        Me.CheckBox4.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox4.Location = New System.Drawing.Point(8, 120)
        Me.CheckBox4.Name = "CheckBox4"
        Me.CheckBox4.TabIndex = 3
        Me.CheckBox4.Text = "PH"
        '
        'CheckBox3
        '
        Me.CheckBox3.Checked = True
        Me.CheckBox3.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox3.Location = New System.Drawing.Point(8, 88)
        Me.CheckBox3.Name = "CheckBox3"
        Me.CheckBox3.TabIndex = 2
        Me.CheckBox3.Text = "Teneur en eau"
        '
        'CheckBox2
        '
        Me.CheckBox2.Checked = True
        Me.CheckBox2.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox2.Location = New System.Drawing.Point(8, 56)
        Me.CheckBox2.Name = "CheckBox2"
        Me.CheckBox2.Size = New System.Drawing.Size(112, 24)
        Me.CheckBox2.TabIndex = 1
        Me.CheckBox2.Text = "Humidité relative"
        '
        'CheckBox1
        '
        Me.CheckBox1.Checked = True
        Me.CheckBox1.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox1.Location = New System.Drawing.Point(8, 24)
        Me.CheckBox1.Name = "CheckBox1"
        Me.CheckBox1.TabIndex = 0
        Me.CheckBox1.Text = "Température"
        '
        'GroupBox3
        '
        Me.GroupBox3.Controls.Add(Me.Label13)
        Me.GroupBox3.Controls.Add(Me.Label12)
        Me.GroupBox3.Controls.Add(Me.Label11)
        Me.GroupBox3.Controls.Add(Me.Label10)
        Me.GroupBox3.Controls.Add(Me.TextBox4)
        Me.GroupBox3.Controls.Add(Me.TextBox3)
        Me.GroupBox3.Controls.Add(Me.Label9)
        Me.GroupBox3.Controls.Add(Me.Label8)
        Me.GroupBox3.Controls.Add(Me.Label7)
        Me.GroupBox3.Controls.Add(Me.Label6)
        Me.GroupBox3.Controls.Add(Me.Label5)
        Me.GroupBox3.Controls.Add(Me.ComboBox1)
        Me.GroupBox3.Location = New System.Drawing.Point(200, 136)
        Me.GroupBox3.Name = "GroupBox3"
        Me.GroupBox3.Size = New System.Drawing.Size(440, 120)
        Me.GroupBox3.TabIndex = 1
        Me.GroupBox3.TabStop = False
        Me.GroupBox3.Text = "Résistance des aciers aux ions chlorures"
        '
        'Label13
        '
        Me.Label13.Location = New System.Drawing.Point(32, 96)
        Me.Label13.Name = "Label13"
        Me.Label13.Size = New System.Drawing.Size(152, 16)
        Me.Label13.TabIndex = 11
        Me.Label13.Text = "Label13"
        '
        'Label12
        '
        Me.Label12.Location = New System.Drawing.Point(32, 72)
        Me.Label12.Name = "Label12"
        Me.Label12.Size = New System.Drawing.Size(152, 16)
        Me.Label12.TabIndex = 10
        Me.Label12.Text = "Label12"
        '
        'Label11
        '
        Me.Label11.Font = New System.Drawing.Font("Symbol", 8.25!, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, CType(2, Byte))
        Me.Label11.Location = New System.Drawing.Point(8, 96)
        Me.Label11.Name = "Label11"
        Me.Label11.Size = New System.Drawing.Size(24, 16)
        Me.Label11.TabIndex = 9
        Me.Label11.Text = "x ="
        '
        'Label10
        '
        Me.Label10.Font = New System.Drawing.Font("Symbol", 8.25!, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, CType(2, Byte))
        Me.Label10.Location = New System.Drawing.Point(8, 72)
        Me.Label10.Name = "Label10"
        Me.Label10.Size = New System.Drawing.Size(24, 16)
        Me.Label10.TabIndex = 8
        Me.Label10.Text = "l ="
        '
        'TextBox4
        '
        Me.TextBox4.Location = New System.Drawing.Point(216, 48)
        Me.TextBox4.Name = "TextBox4"
        Me.TextBox4.TabIndex = 7
        Me.TextBox4.Text = ""
        '
        'TextBox3
        '
        Me.TextBox3.Location = New System.Drawing.Point(216, 24)
        Me.TextBox3.Name = "TextBox3"
        Me.TextBox3.TabIndex = 6
        Me.TextBox3.Text = ""
        '
        'Label9
        '
        Me.Label9.Location = New System.Drawing.Point(320, 48)
        Me.Label9.Name = "Label9"
        Me.Label9.Size = New System.Drawing.Size(104, 16)
        Me.Label9.TabIndex = 5
        Me.Label9.Text = "% masse de ciment"
        '
        'Label8
        '
        Me.Label8.Location = New System.Drawing.Point(320, 24)
        Me.Label8.Name = "Label8"
        Me.Label8.Size = New System.Drawing.Size(104, 16)
        Me.Label8.TabIndex = 4
        Me.Label8.Text = "% masse de ciment"
        '
        'Label7
        '
        Me.Label7.Location = New System.Drawing.Point(8, 48)
        Me.Label7.Name = "Label7"
        Me.Label7.Size = New System.Drawing.Size(224, 16)
        Me.Label7.TabIndex = 3
        Me.Label7.Text = "teneur en ions chlorures libres (écart-type)"
        '
        'Label6
        '
        Me.Label6.Location = New System.Drawing.Point(8, 24)
        Me.Label6.Name = "Label6"
        Me.Label6.Size = New System.Drawing.Size(216, 16)
        Me.Label6.TabIndex = 2
        Me.Label6.Text = "teneur en ions chlorures libres (moyenne)"
        '
        'Label5
        '
        Me.Label5.Location = New System.Drawing.Point(200, 96)
        Me.Label5.Name = "Label5"
        Me.Label5.Size = New System.Drawing.Size(96, 16)
        Me.Label5.TabIndex = 1
        Me.Label5.Text = "Catégorie d'acier"
        '
        'ComboBox1
        '
        Me.ComboBox1.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList
        Me.ComboBox1.Items.AddRange(New Object() {"acier usuel", "acier Top12", "acier inoxydable", "autres"})
        Me.ComboBox1.Location = New System.Drawing.Point(312, 88)
        Me.ComboBox1.Name = "ComboBox1"
        Me.ComboBox1.Size = New System.Drawing.Size(121, 21)
        Me.ComboBox1.TabIndex = 0
        '
        'GroupBox2
        '
        Me.GroupBox2.Controls.Add(Me.TextBox2)
        Me.GroupBox2.Controls.Add(Me.TextBox1)
        Me.GroupBox2.Controls.Add(Me.Label4)
        Me.GroupBox2.Controls.Add(Me.Label3)
        Me.GroupBox2.Controls.Add(Me.Label2)
        Me.GroupBox2.Controls.Add(Me.Label1)
        Me.GroupBox2.Controls.Add(Me.CheckBox8)
        Me.GroupBox2.Location = New System.Drawing.Point(200, 8)
        Me.GroupBox2.Name = "GroupBox2"
        Me.GroupBox2.Size = New System.Drawing.Size(440, 112)
        Me.GroupBox2.TabIndex = 2
        Me.GroupBox2.TabStop = False
        Me.GroupBox2.Text = "Enrobage des aciers"
        '
        'TextBox2
        '
        Me.TextBox2.Location = New System.Drawing.Point(192, 56)
        Me.TextBox2.Name = "TextBox2"
        Me.TextBox2.TabIndex = 6
        Me.TextBox2.Text = "23"
        '
        'TextBox1
        '
        Me.TextBox1.Location = New System.Drawing.Point(192, 32)
        Me.TextBox1.Name = "TextBox1"
        Me.TextBox1.TabIndex = 5
        Me.TextBox1.Text = "35"
        '
        'Label4
        '
        Me.Label4.Location = New System.Drawing.Point(312, 56)
        Me.Label4.Name = "Label4"
        Me.Label4.Size = New System.Drawing.Size(24, 16)
        Me.Label4.TabIndex = 4
        Me.Label4.Text = "mm"
        '
        'Label3
        '
        Me.Label3.Location = New System.Drawing.Point(312, 32)
        Me.Label3.Name = "Label3"
        Me.Label3.Size = New System.Drawing.Size(24, 23)
        Me.Label3.TabIndex = 3
        Me.Label3.Text = "mm"
        '
        'Label2
        '
        Me.Label2.Location = New System.Drawing.Point(8, 56)
        Me.Label2.Name = "Label2"
        Me.Label2.Size = New System.Drawing.Size(152, 16)
        Me.Label2.TabIndex = 2
        Me.Label2.Text = "position de l'acier (écart-type)"
        '
        'Label1
        '
        Me.Label1.Location = New System.Drawing.Point(8, 32)
        Me.Label1.Name = "Label1"
        Me.Label1.Size = New System.Drawing.Size(152, 16)
        Me.Label1.TabIndex = 1
        Me.Label1.Text = "position de l'acier (moyenne)"
        '
        'CheckBox8
        '
        Me.CheckBox8.Checked = True
        Me.CheckBox8.CheckState = System.Windows.Forms.CheckState.Checked
        Me.CheckBox8.Location = New System.Drawing.Point(16, 80)
        Me.CheckBox8.Name = "CheckBox8"
        Me.CheckBox8.Size = New System.Drawing.Size(144, 24)
        Me.CheckBox8.TabIndex = 0
        Me.CheckBox8.Text = "Traitement probabiliste"
        '
        'Button1
        '
        Me.Button1.Location = New System.Drawing.Point(8, 264)
        Me.Button1.Name = "Button1"
        Me.Button1.TabIndex = 3
        Me.Button1.Text = "Traitement"
        '
        'Button2
        '
        Me.Button2.Location = New System.Drawing.Point(200, 264)
        Me.Button2.Name = "Button2"
        Me.Button2.TabIndex = 4
        Me.Button2.Text = "Analyse"
        '
        'Button3
        '
        Me.Button3.Location = New System.Drawing.Point(568, 264)
        Me.Button3.Name = "Button3"
        Me.Button3.TabIndex = 5
        Me.Button3.Text = "Sortir"
        '
        'frmProb
        '
        Me.AutoScaleBaseSize = New System.Drawing.Size(5, 13)
        Me.AutoScroll = True
        Me.ClientSize = New System.Drawing.Size(648, 296)
        Me.Controls.Add(Me.Button3)
        Me.Controls.Add(Me.Button2)
        Me.Controls.Add(Me.Button1)
        Me.Controls.Add(Me.GroupBox2)
        Me.Controls.Add(Me.GroupBox3)
        Me.Controls.Add(Me.GroupBox1)
        Me.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle
        Me.Name = "frmProb"
        Me.Text = "Approche probabiliste"
        Me.GroupBox1.ResumeLayout(False)
        Me.GroupBox3.ResumeLayout(False)
        Me.GroupBox2.ResumeLayout(False)
        Me.ResumeLayout(False)

    End Sub

#End Region

    'chargement de la fenêtre
    Private Sub Form1_Load(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MyBase.Load
        ComboBox1.SelectedIndex = 3
        Button2.Visible = False
    End Sub

    'click sur traitement
    Private Sub Button1_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button1.Click
        Dim Data As Object
        Dim La(1, 1, 1, 1, 1, 1) As Double
        Dim Cpd(2, 2, 2, 2) As Double
        Dim Filtre As String
        Dim Index As Short
        Dim Directoire As Boolean
        Dim Titre As String
        Dim Canc As Boolean = False
        Dim OutF As String
        Dim Poub As String
        Dim Test01 As String
        Dim TProb As String
        Dim prefFile As String
        Dim Tprob1 As Short
        Dim Tprob2 As Short
        Dim Tprob3 As Short
        Dim Tprob4 As Short
        Dim Ncol As Short
        Dim i As Short
        Dim j As Short
        Dim k As Short
        Dim l As Short
        Dim m As Short
        Dim n As Short
        Dim o As Short
        Dim p As Short
        Dim Dim1 As String
        Dim Dim2 As String
        Dim iPos As Short
        Dim nFic As Short
        Dim pd As Double
        Dim pd_cont As Double = 0
        Dim cpt1 As Short = 0
        Dim cpt2 As Short = 0

        ''''''''''''''''''''''''''''''''''''''''''''
        Filtre = "txt files (R_*.txt)|R_*.txt|All files (*.*)|*.*"
        Index = CShort(1)
        Directoire = True
        Titre = "Sélectionner un fichier résultat à traiter"
        OpenDialog(OutFile, Canc, Filtre, Index, Directoire, Titre)
        If Canc = True Then GoTo b
        '''''''''''''''''''''''''''''''''''''''''''''
        prefFile = OutFile
        FileOnly(OutFile)
        prefFile = Microsoft.VisualBasic.Left(prefFile, CInt(Len(prefFile) - Len(OutFile)))
        OutFile = Microsoft.VisualBasic.Left(OutFile, Len(OutFile) - 4)
        TProb = Microsoft.VisualBasic.Right(OutFile, 4)
        OutFile = Microsoft.VisualBasic.Left(OutFile, Len(OutFile) - 4)
        OutFile = Microsoft.VisualBasic.Right(OutFile, Len(OutFile) - 2)
        Dim1 = "_"
        iPos = InStr(1, OutFile, Dim1, CompareMethod.Text) 'position de "_"
        OutFile = Microsoft.VisualBasic.Right(OutFile, Len(OutFile) - iPos)
        Tprob1 = CShort(Microsoft.VisualBasic.Left(TProb, 1))
        Tprob4 = CShort(Microsoft.VisualBasic.Right(TProb, 1))
        Tprob2 = (TProb - 1000 * Tprob1 - tprob4) / 10
        Tprob3 = CShort(Microsoft.VisualBasic.Right(Tprob2, 1))
        Tprob2 = CShort(Microsoft.VisualBasic.Left(Tprob2, 1))
        If Tprob1 = 0 Then
            Tprob1 = 1
        Else
            Tprob1 = 2
            cpt1 = 2
        End If
        If Tprob2 = 0 Then
            Tprob2 = 1
        Else
            Tprob2 = 2
            cpt1 = cpt1 * 2
            If cpt1 = 0 Then cpt1 = 2
        End If
        If Tprob3 = 0 Then
            Tprob3 = 1
        Else
            Tprob3 = 2
            cpt1 = cpt1 * 2
            If cpt1 = 0 Then cpt1 = 2
        End If
        If Tprob4 = 0 Then
            Tprob4 = 1
        Else
            Tprob4 = 2
            cpt1 = cpt1 * 2
            If cpt1 = 0 Then cpt1 = 2
        End If

        Files = 0
        If CheckBox1.Checked = True Then
            Files = Files + 1
            PreFile(Files) = "R_T_"
        End If
        If CheckBox2.Checked = True Then
            Files = Files + 1
            PreFile(Files) = "R_H_"
        End If
        If CheckBox3.Checked = True Then
            Files = Files + 1
            PreFile(Files) = "R_W_"
        End If
        If CheckBox4.Checked = True Then
            Files = Files + 1
            PreFile(Files) = "R_PH_"
        End If
        If CheckBox6.Checked = True Then
            Files = Files + 1
            PreFile(Files) = "R_CL_"
        End If
        If CheckBox7.Checked = True Then
            Files = Files + 1
            PreFile(Files) = "R_CT_"
        End If
        If CheckBox5.Checked = True Then
            Files = Files + 1
            PreFile(Files) = "R_Carb_"
        End If

        For o = 1 To Files
            pd_cont = 0
            cpt2 = 0
            For n = 1 To Tprob4
                For i = 1 To Tprob3
                    For j = 1 To Tprob2
                        For k = 1 To Tprob1

                            cpt2 = cpt2 + 1
                            If Tprob1 = 1 Then k = 0
                            If Tprob2 = 1 Then j = 0
                            If Tprob3 = 1 Then i = 0
                            If Tprob4 = 1 Then n = 0
                            OutF = prefFile & PreFile(o) & OutFile & k & j & i & n & ".txt"
                            If k = 0 Then k = k + 1
                            If j = 0 Then j = j + 1
                            If i = 0 Then i = i + 1
                            If n = 0 Then n = n + 1
                            nFic = FreeFile()
                            FileOpen(nFic, OutF, OpenMode.Input, OpenAccess.Read, OpenShare.Shared)

                            Input(nFic, Dim1) 'nbre de ligne
                            Input(nFic, Dim2) 'nbre de colonnes
                            Dim1 = "_"
                            iPos = InStr(1, Dim2, Dim1, CompareMethod.Text) 'position de "_"
                            If cpt2 = 1 Then Nline = CShort(Mid(Dim2, 1, iPos - 1))
                            Dim2 = CStr(Mid(Dim2, iPos + 1))
                            Dim1 = " "
                            iPos = 1
                            Do While iPos = 1
                                iPos = InStr(1, Dim2, Dim1, CompareMethod.Text)
                                Dim2 = CStr(Mid(Dim2, iPos + 1))
                                iPos = InStr(1, Dim2, Dim1, CompareMethod.Text)
                            Loop
                            If cpt2 = 1 Then Ncol = CShort(Mid(Dim2, 1, iPos + 1))
                            If cpt2 = 1 And o = 1 Then   'lecture de la première ligne
                                ReDim Labe(Ncol - 2)
                                ReDim Lambda(Nline + 2, Ncol - 2, Files)
                                ReDim Ksi(Nline + 2, Ncol - 2, Files)
                                ReDim La(Nline + 2, Ncol - 2, 2, 2, 2, 2)
                                Nc = Ncol - 2
                                If PreFile(o) = "R_H_" Or PreFile(o) = "R_W_" Then
                                    ReDim Lambda(Nline + 2, Ncol - 3, Files)
                                    ReDim Ksi(Nline + 2, Ncol - 3, Files)
                                    ReDim La(Nline + 2, Ncol - 3, 2, 2, 2, 2)
                                    ReDim Labe(Ncol - 3)
                                    Nc = Ncol - 3
                                End If
                                If PreFile(o) = "R_Carb_" Then
                                    ReDim Lambda(Nline + 2, 4, Files)
                                    ReDim Ksi(Nline + 2, 4, Files)
                                    ReDim Labe(4)
                                    Nc = 4
                                End If
                                If PreFile(o) = "R_T_" Or PreFile(o) = "R_PH_" Or PreFile(o) = "R_Carb_" Or PreFile(o) = "R_CL_" Or PreFile(o) = "R_CT_" Then p = 2
                                If PreFile(o) = "R_H_" Or PreFile(o) = "R_W_" Then p = 3
                                If PreFile(o) = "R_Carb_" Then
                                    Ncol = 5
                                    p = 2
                                End If
                                For l = 2 To Ncol - 1
                                    Input(nFic, Data)
                                    If l <= Ncol - p Then
                                        Length = CSng(Data)
                                        Labe(l) = Data
                                    End If
                                Next l
                            Else
                                If PreFile(o) = "R_Carb_" Then Ncol = 5
                                For l = 2 To Ncol - 1
                                    Input(nFic, Poub)
                                Next l
                            End If

                            Dim1 = "_"  'poids des fonctions probabilistes
                            iPos = InStr(1, Dim2, Dim1, CompareMethod.Text) 'position de "_"
                            Dim2 = Microsoft.VisualBasic.Right(Dim2, Len(Dim2) - iPos)
                            Dim1 = " "
                            iPos = 1
                            Do While iPos = 1
                                iPos = InStr(1, Dim2, Dim1, CompareMethod.Text)
                                Dim2 = CStr(Mid(Dim2, iPos + 1))
                                iPos = InStr(1, Dim2, Dim1, CompareMethod.Text)
                            Loop
                            pd = CDbl(Dim2)
                            pd_cont = pd_cont + pd  'contrôle
                            Cpd(k, j, i, n) = pd
                            If cpt2 = cpt1 And o = 1 Then Test01 = "La somme des poids des fonctions probabilistes vaut :" & pd_cont
                            If cpt2 = cpt1 And CSng(pd_cont) = 1 And o = 1 Then Test01 = "Réussite du contrôle la somme des poids des fonctions probabilistes"
                            If cpt2 = cpt1 And o <> 1 Then
                                If CSng(pd_cont) = 1 Then
                                    Test01 = Test01 & " ok"
                                Else
                                    Test01 = Test01 & " " & pd_cont
                                End If
                            End If

                            If PreFile(o) = "R_T_" Or PreFile(o) = "R_PH_" Or PreFile(o) = "R_Carb_" Or PreFile(o) = "R_CL_" Or PreFile(o) = "R_CT_" Then p = 2
                            If PreFile(o) = "R_H_" Or PreFile(o) = "R_W_" Then p = 3
                            If PreFile(o) = "R_Carb_" Then
                                Ncol = 5
                                p = 2
                            End If
                            For l = 1 To Nline + 2         'lecture du fichier résultat
                                For m = 0 To Ncol - 1
                                    Try
                                        Input(nFic, Data)
                                    Catch ex As Exception When Data = ""
                                        Exit For
                                        Exit For
                                    End Try
                                    If m <= Ncol - p Then
                                        If m <= 1 Then
                                            Lambda(l, m, o) = CDbl(Data)
                                        Else
                                            If Data <= 0 Then Data = 1.0E-300 'annuler les problèmes numériques
                                            Lambda(l, m, o) = Lambda(l, m, o) + pd * System.Math.Log(CDbl(Data))
                                            La(l, m, k, j, i, n) = CDbl(Data)
                                        End If
                                    End If
                                Next m
                            Next l
                            FileClose(nFic)

                            If cpt2 = cpt1 And o = Files Then
                                MsgBox(Test01, MsgBoxStyle.Information And MsgBoxStyle.OKOnly, "Avertissement")
                                Button2.Visible = True
                                GroupBox1.Enabled = False
                                Button1.Enabled = False
                            End If
                        Next k
                    Next j
                Next i
            Next n
            For n = 1 To Tprob4
                For i = 1 To Tprob3
                    For j = 1 To Tprob2
                        For k = 1 To Tprob1
                            For l = 1 To Nline + 2
                                For m = 2 To Ncol - p
                                    Ksi(l, m, o) = Ksi(l, m, o) + Cpd(k, j, i, n) * (System.Math.Log(CDbl(La(l, m, k, j, i, n))) - Lambda(l, m, o)) ^ 2
                                Next m
                            Next l
                        Next k
                    Next j
                Next i
            Next n
        Next o
b:
    End Sub

    'click sur sortir
    Private Sub Button3_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button3.Click
        GroupBox1.Enabled = True
        Button1.Enabled = True
        Me.Hide()
    End Sub

    'changement du type d'acier
    Private Sub ComboBox1_SelectedIndexChanged(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles ComboBox1.SelectedIndexChanged
        Dim a As Decimal
        Dim b As Decimal
        Select Case ComboBox1.Text
            Case "acier usuel"
                TextBox3.Text = CStr(0.536239738583839)
                TextBox4.Text = CStr(0.761717254531377)
                TextBox3.Enabled = False
                TextBox4.Enabled = False
            Case "acier Top12"
                TextBox3.Text = CStr(1.3405993464596)
                TextBox4.Text = CStr(0.761717254531377)
                TextBox3.Enabled = False
                TextBox4.Enabled = False
            Case "acier inoxydable"
                TextBox3.Text = CStr(2.68119869291919)
                TextBox4.Text = CStr(0.761717254531377)
                TextBox3.Enabled = False
                TextBox4.Enabled = False
            Case "autres"
                TextBox3.Text = ""
                TextBox4.Text = ""
                Label12.Text = ""
                Label13.Text = ""
                TextBox3.Enabled = True
                TextBox4.Enabled = True
        End Select
    End Sub

    'changement de la moyenne et écart-type dans le type d'acier
    Private Sub AcierChange(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles TextBox3.TextChanged, TextBox4.TextChanged
        If TextBox3.Text = "" Or TextBox4.Text = "" Then Exit Sub
        Label12.Text = System.Math.Log(CDbl(TextBox3.Text) ^ 2 * (CDbl(TextBox3.Text) ^ 2 + CDbl(TextBox4.Text) ^ 2) ^ 0.5)
        Label13.Text = (System.Math.Log(CDbl(TextBox4.Text) ^ 2 / CDbl(TextBox3.Text) ^ 2 + 1)) ^ 0.5
    End Sub

    'click sur analyse
    Private Sub Button2_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button2.Click
        Dim LAen(Nline + 2, Files) As Double
        Dim KSen(Nline + 2, Files) As Double
        Dim Pf(Nline + 2) As Double
        Dim Userinput2 As String
        Dim PlFile As String
        Dim OT As String
        Dim i As Short
        Dim j As Short
        Dim k As Short
        Dim o As Short
        Dim Men As Single
        Dim ETen As Single
        Dim PosX As Single = 0
        Dim PosmX As Single = 0
        Dim PosXi As Single
        Dim deltaX1 As Single = 0.01
        Dim deltaX2 As Single = 0.1
        Dim deltaX3 As Single = 0.01
        Dim deltaCL2 As Single = 0.005
        Dim Ca As Boolean = False
        Dim Var1 As Double
        Dim Var2 As Double
        Dim Var3 As Double
        Dim Var4 As Double
        Dim Var5 As Double
        Dim Var6 As Double
        Dim Var7 As Double
        Dim Var8 As Double
        Dim Var9 As Double
        Dim nfic As Short
        Dim Spos As Boolean
        Dim Spos_old As Boolean
        Dim Request As Short
        Dim cpt As Short
        Dim val As Short

        Dim Filtre As String
        Dim Index As Short
        Dim Canc As Boolean = False
        Dim Directoire As Boolean
        Dim Titre As String

        If PreFile(1) <> "R_Carb_" Then
            Request = MsgBox("A partir d'un fichier existant LP_ ?", MsgBoxStyle.YesNo, "Question")
            If Request = MsgBoxResult.No Then
                If PreFile(Files) = "R_Carb_" Then Ca = True
                If Ca = True Then Files = Files - 1 'enlever la partie carbonatation
                If CheckBox8.Checked = True Then
                    PosX = 0
                    PosmX = 0
                    Men = CSng(TextBox1.Text)
                    ETen = CSng(TextBox2.Text)
                    Var6 = 0
                    Do While PosX < Length

                        PosX = PosX + deltaX1
                        If PosX >= Length Then PosX = Length
                        deltaX1 = PosX - PosmX
                        Var1 = System.Math.Exp(-0.5 * ((PosmX - Men) / ETen) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * ETen)
                        Var2 = System.Math.Exp(-0.5 * ((PosX - Men) / ETen) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * ETen)
                        Var3 = (Var1 + Var2) * deltaX1 / 2          'Surface par la méthode des trapèzes
                        Var6 = Var3 + Var6

                        For i = 2 To Nc
                            PosXi = (PosX + PosmX) / 2
                            If Labe(i) > PosXi Then Exit For
                        Next i
                        For o = 1 To Files
                            For j = 1 To Nline + 2
                                Var4 = (Lambda(j, i, o) - Lambda(j, i - 1, o)) * (PosXi - Labe(i)) / (Labe(i) - Labe(i - 1)) + Lambda(j, i, o)
                                Var5 = (Ksi(j, i, o) - Ksi(j, i - 1, o)) * (PosXi - Labe(i)) / (Labe(i) - Labe(i - 1)) + Ksi(j, i, o)
                                LAen(j, o) = LAen(j, o) + Var4 * Var3
                                KSen(j, o) = KSen(j, o) + Var5 * Var3
                                If PosX = Length Then
                                    LAen(j, o) = LAen(j, o) / Var6
                                    KSen(j, o) = KSen(j, o) / Var6
                                End If
                            Next j
                        Next o
                        PosmX = PosX
                    Loop
                Else
                    Userinput2 = InputBox("Approche déterministe de l'enrobage, profondeur d'enrobage en [mm], (0," & Length & ") : ", "Profondeur d'enrobage des aciers d'armature", "35")
                    Msg_noEntry(Userinput2, Canc)
                    If Canc = False Then Msg_noNumeric(Userinput2, Canc)
                    If Userinput2 < 0 Or Userinput2 > Length Then Canc = True
                    If Canc = True Then Exit Sub
                    PosXi = CShort(Userinput2)

                    For i = 2 To Nc
                        If Labe(i) > PosXi Then Exit For
                    Next i
                    For o = 1 To Files
                        For j = 1 To Nline + 2
                            LAen(j, o) = (Lambda(j, i, o) - Lambda(j, i - 1, o)) * (PosXi - Labe(i)) / (Labe(i) - Labe(i - 1)) + Lambda(j, i, o)
                            KSen(j, o) = (Ksi(j, i, o) - Ksi(j, i - 1, o)) * (PosXi - Labe(i)) / (Labe(i) - Labe(i - 1)) + Ksi(j, i, o)
                        Next j
                    Next o
                End If

                OT = "LP_" & OutFile        'enregistrement intermédiaire
                ''''''''''''''''''''''''''''''''''
                Filtre = "txt files (LP_*.txt)|LP_*.txt|All files (*.*)|*.*"
                Index = CShort(1)
                Directoire = True
                Titre = "Paramètres des lois probabilistes"
                SaveDialog(OT, Canc, Filtre, Index, Directoire, Titre)
                If Canc = True Then GoTo b
                ''''''''''''''''''''''''''''''''''
                nfic = CShort(FreeFile())
                FileOpen(CInt(nfic), OT, OpenMode.Output)
                If CheckBox8.Checked = True Then
                    PrintLine(CInt(nfic), "Approche probabiliste , enrobage moyen ", Men, " écart-type de l'enrobage ", ETen)
                Else
                    PrintLine(CInt(nfic), "Approche déterministe de l'enrobage , valeur ", PosXi)
                End If
                If Ca = True Then
                    Print(nfic, "temps, temps,")
                    For o = 1 To Files + 1
                        Print(nfic, PreFile(o), ",", PreFile(o), ",")
                    Next o
                    PrintLine(nfic, PreFile(Files + 1), ",", PreFile(Files + 1), ",")
                    Print(nfic, "années, jours,")
                    For o = 1 To Files + 2
                        Print(nfic, "lambda,ksi,")
                    Next
                    PrintLine(nfic, "")
                    For j = 1 To Nline + 2
                        Print(nfic, Lambda(j, 0, 1), ",", Lambda(j, 1, 1), ",")
                        For o = 1 To Files
                            Print(nfic, LAen(j, o), ",", KSen(j, o), ",")
                        Next o
                        PrintLine(nfic, Lambda(j, 2, Files + 1), ",", Ksi(j, 2, Files + 1), ",", Lambda(j, 3, Files + 1), ",", Ksi(j, 3, Files + 1), ",")
                    Next j
                Else
                    Print(nfic, "temps, temps,")
                    For o = 1 To Files
                        Print(nfic, PreFile(o), ",", PreFile(o), ",")
                    Next o
                    PrintLine(nfic, "")
                    Print(nfic, "années, jours,")
                    For o = 1 To Files
                        Print(nfic, "lambda,ksi,")
                    Next
                    PrintLine(nfic, "")
                    For j = 1 To Nline + 2
                        Print(nfic, Lambda(j, 0, 1), ",", Lambda(j, 1, 1), ",")
                        For o = 1 To Files
                            Print(nfic, LAen(j, o), ",", KSen(j, o), ",")
                        Next o
                        PrintLine(nfic, "")
                    Next j
                End If
                FileClose(nfic)
            Else
                ''''''''''''''''''''''''''''''''''''''''''''
                Filtre = "txt files (LP_*.txt)|LP_*.txt|All files (*.*)|*.*"
                Index = CShort(1)
                Directoire = True
                Titre = "Sélectionner un fichier résultat à traiter"
                OpenDialog(PlFile, Canc, Filtre, Index, Directoire, Titre)
                If Canc = True Then GoTo b
                '''''''''''''''''''''''''''''''''''''''''''''
                nfic = FreeFile()
                FileOpen(nfic, PlFile, OpenMode.Input, OpenAccess.Read, OpenShare.Shared)
                Input(nfic, Userinput2)
                Input(nfic, Userinput2)
                Input(nfic, Userinput2)
                Input(nfic, Userinput2)
                Input(nfic, Userinput2)
                cpt = 1
                val = 0
                Do While Microsoft.VisualBasic.Left(Userinput2, 2) = "R_"
                    Input(nfic, Userinput2)
                    If Userinput2 = "R_Carb_" Then
                        val = 2
                        Ca = True
                    End If
                    cpt = cpt + 1
                Loop
                cpt = cpt - 1
                Files = cpt / 2 - val
                If val = 0 Then Ca = False
                If Ca = True Then
                    For o = 1 To (Files + 3) * 2 + 1
                        Input(nfic, Userinput2)
                    Next o
                    j = 0
                    Do While Userinput2 <> "-10"
                        j = j + 1
                        Try
                            Input(nfic, Userinput2)
                        Catch ex As Exception When Userinput2 = ""
                            Exit Do
                        End Try
                        Input(nfic, Userinput2)
                        For o = 1 To Files
                            Input(nfic, LAen(j, o))
                            Input(nfic, KSen(j, o))
                        Next o
                        Input(nfic, Lambda(j, 2, Files + 1))
                        Input(nfic, Ksi(j, 2, Files + 1))
                        Input(nfic, Lambda(j, 3, Files + 1))
                        Input(nfic, Ksi(j, 3, Files + 1))
                        Input(nfic, Userinput2)
                    Loop
                Else
                    For o = 1 To Files * 2 + 1
                        Input(nfic, Userinput2)
                    Next o
                    j = 0
                    Do While Userinput2 <> "-10"
                        j = j + 1
                        Try
                            Input(nfic, Userinput2)
                        Catch ex As Exception When Userinput2 = ""
                            Exit Do
                        End Try
                        Input(nfic, Userinput2)
                        For o = 1 To Files
                            Input(nfic, LAen(j, o))
                            Input(nfic, KSen(j, o))
                        Next o
                        Input(nfic, Userinput2)
                    Loop
                End If
                FileClose(nfic)
            End If

b:
            If TextBox3.Text <> "" And TextBox4.Text <> "" Then         'calcul de la probabilité d'initiation de la corrosion due à la présence de ions chlorures
                OT = "PFCL_" & OutFile        'enregistrement des probabilité d'initiation de la corrosion due à la présence des ions chlorures
                ''''''''''''''''''''''''''''''''''
                Filtre = "txt files (PFCL_*.txt)|PFCL_*.txt|All files (*.*)|*.*"
                Index = CShort(1)
                Directoire = True
                Titre = "Probabilités d'initiation de la corrosion (Cl-)"
                SaveDialog(OT, Canc, Filtre, Index, Directoire, Titre)
                If Canc = True Then GoTo d
                ''''''''''''''''''''''''''''''''''
                For o = 1 To Files
                    If PreFile(o) = "R_CL_" Then
                        PosX = 0
                        Exit For
                    Else
                        PosX = 1
                    End If
                Next o
                If PosX = 0 Then
                    For j = 1 To Nline + 2
                        PosX = 0
                        Var7 = 1
                        Var5 = 1.0E-300
                        Var8 = 0
                        Do While PosX < 2 Or Var7 >= 0.00001
                            PosX = PosX + deltaCL2
                            Var1 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX) - LAen(j, o)) / KSen(j, o)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * KSen(j, o) * PosX)            'fs
                            Var3 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX) - CDbl(Label12.Text)) / CDbl(Label13.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(Label13.Text) * PosX)
                            Var6 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX - deltaCL2 / 2) - CDbl(Label12.Text)) / CDbl(Label13.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(Label13.Text) * (PosX - deltaCL2 / 2))
                            Var2 = System.Math.Exp(-0.5 * ((System.Math.Log(Var5) - CDbl(Label12.Text)) / CDbl(Label13.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(Label13.Text) * Var5)
                            Var2 = (Var2 + Var6) * deltaCL2 / 4         'méthode des trapèzes
                            Var3 = (Var3 + Var6) * deltaCL2 / 4 + Var2
                            Var9 = Var1 * (Var8 + Var3) 'partie de Pf sans l'intégration de xi+1

                            Var1 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX - deltaCL2 / 2) - LAen(j, o)) / KSen(j, o)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * KSen(j, o) * (PosX - deltaCL2 / 2))    'fs
                            Var6 = Var1 * (Var8 + Var2) 'partie de Pf sans l'intégration de (xi+1 + xi)/2
                            If PosX <> deltaCL2 Then        'formule de simpson
                                Var7 = (PosX - Var5) * (Var9 + 4 * Var6 + Var4) / 6
                                Pf(j) = Pf(j) + Var7
                            End If
                            Var8 = Var8 + Var3
                            Var4 = Var9
                            Var5 = PosX
                        Loop
                    Next j

                    nfic = CShort(FreeFile())
                    FileOpen(CInt(nfic), OT, OpenMode.Output)
                    PrintLine(CInt(nfic), "Probabilité d'initiation de la corrosion due à la présence de ions chlorures")
                    PrintLine(nfic, "temps, temps, Pf,")
                    PrintLine(nfic, "années, jours,")
                    For j = 1 To Nline + 2
                        PrintLine(nfic, Lambda(j, 0, 1), ",", Lambda(j, 1, 1), ",", Pf(j), ",")
                    Next j
                    FileClose(nfic)

                    OT = "BCL_" & OutFile        'enregistrement de l'indice de fiabilité pour la corrosion par les chlorures
                    ''''''''''''''''''''''''''''''''''
                    Filtre = "txt files (BCL_*.txt)|BCL_*.txt|All files (*.*)|*.*"
                    Index = CShort(1)
                    Directoire = True
                    Titre = "Indice de fiabilité (Cl-)"
                    SaveDialog(OT, Canc, Filtre, Index, Directoire, Titre)
                    If Canc = True Then GoTo d
                    ''''''''''''''''''''''''''''''''''
                    For j = 1 To Nline + 2  'calcul de l'indice de fiabilité bêta
                        Var1 = -4           'beta
                        Var2 = 4            'delta beta
                        Var4 = 10
                        Spos = False
                        Spos_old = False
                        Do While System.Math.Abs(Var2) > 0.001
                            If Var4 <> 10 Then Var1 = Var1 + Var2
                            If Var1 >= 0 Then
                                Var3 = 1 / (1 + 0.3275911 * Var1)
                                Var3 = 1 - (0.254829592 * Var3 - 0.284496736 * Var3 ^ 2 + 1.421413741 * Var3 ^ 3 - 1.453152027 * Var3 ^ 4 + 1.061405429 * Var3 ^ 5) * System.Math.Exp(-Var1 ^ 2)
                                Var3 = 0.5 * (1 - Var3)         'pf approximé
                            Else
                                Var3 = 1 / (1 + 0.3275911 * (-Var1))
                                Var3 = 1 - (0.254829592 * Var3 - 0.284496736 * Var3 ^ 2 + 1.421413741 * Var3 ^ 3 - 1.453152027 * Var3 ^ 4 + 1.061405429 * Var3 ^ 5) * System.Math.Exp(-Var1 ^ 2)
                                Var3 = 0.5 * (1 + Var3)         'pf approximé
                            End If
                            Var4 = Pf(j) - Var3
                            If Var4 < 0 Then
                                Spos = False
                            Else
                                Spos = True
                            End If
                            If Spos <> Spos_old Then
                                Var2 = -Var2 / 10
                            End If
                            Spos_old = Spos
                            If Pf(j) = 0 Then
                                Var1 = 20
                                Exit Do
                            End If
                            If Var1 < -20 Then
                                Var1 = -20
                                Exit Do
                            End If
                        Loop
                        KSen(j, 0) = Var1           'indice de fiabilité bêta
                    Next j

                    nfic = CShort(FreeFile())
                    FileOpen(CInt(nfic), OT, OpenMode.Output)
                    PrintLine(CInt(nfic), "Indice de fiabilité pour la corrosion due à la présence de ions chlorures")
                    PrintLine(nfic, "temps, temps, Bêta,")
                    PrintLine(nfic, "années, jours,")
                    For j = 1 To Nline + 2
                        PrintLine(nfic, Lambda(j, 0, 1), ",", Lambda(j, 1, 1), ",", KSen(j, 0), ",")
                        LAen(j, 0) = 0
                        LAen(j, 1) = 0
                    Next j
                    FileClose(nfic)
d:
                End If
            End If
        Else
                Ca = True
            End If

            If Ca = True And CheckBox8.Checked = True Then          'approche probabiliste pour la carbonatation
            OT = "PFCa_" & OutFile        'enregistrement des probabilité d'initiation de la corrosion due à la carbonatation
            ''''''''''''''''''''''''''''''''''
            Filtre = "txt files (PFCa_*.txt)|PFCa_*.txt|All files (*.*)|*.*"
            Index = CShort(1)
            Directoire = True
            Titre = "Probabilités d'initiation de la corrosion (carb)"
            SaveDialog(OT, Canc, Filtre, Index, Directoire, Titre)
            If Canc = True Then GoTo f
            ''''''''''''''''''''''''''''''''''
            If PreFile(1) = "R_Carb_" Then Files = Files - 1
            For j = 1 To Nline + 2
                PosX = 0       'carbonatation depuis la gauche
                Var7 = 1
                Var5 = 0
                Var8 = 0
                If j = 1 Then PosmX = deltaX2
                deltaX2 = PosmX
                Do While PosX < Length
                    PosX = PosX + deltaX2
                    If PosX > Length Then
                        deltaX2 = Length - PosX
                        PosX = Length
                    End If
                    Var1 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX) - Lambda(j, 2, Files + 1)) / Ksi(j, 2, Files + 1)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * Ksi(j, 2, Files + 1) * PosX)            'fs
                    Var3 = System.Math.Exp(-0.5 * ((PosX - CDbl(TextBox1.Text)) / CDbl(TextBox2.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(TextBox2.Text))
                    Var6 = System.Math.Exp(-0.5 * (((PosX - deltaX2 / 2) - CDbl(TextBox1.Text)) / CDbl(TextBox2.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(TextBox2.Text))
                    Var2 = System.Math.Exp(-0.5 * ((Var5 - CDbl(TextBox1.Text)) / CDbl(TextBox2.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(TextBox2.Text))
                    Var2 = (Var2 + Var6) * deltaX2 / 4         'méthode des trapèzes
                    Var3 = (Var3 + Var6) * deltaX2 / 4 + Var2
                    Var9 = Var1 * (Var8 + Var3) 'partie de Pf sans l'intégration de xi+1

                    Var1 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX - deltaX2 / 2) - Lambda(j, 2, Files + 1)) / Ksi(j, 2, Files + 1)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * Ksi(j, 2, Files + 1) * (PosX - deltaCL2 / 2))    'fs
                    Var6 = Var1 * (Var8 + Var2) 'partie de Pf sans l'intégration de (xi+1 + xi)/2
                    If PosX <> deltaX2 Then        'formule de simpson
                        Var7 = (PosX - Var5) * (Var9 + 4 * Var6 + Var4) / 6
                        LAen(j, 0) = LAen(j, 0) + Var7
                    End If
                    Var8 = Var8 + Var3
                    Var4 = Var9
                    Var5 = PosX
                Loop
                LAen(j, 0) = LAen(j, 0) / Var8  'correction due au troncage de l'épaisseur du mur

                'PosX = Length       'carbonatation depuis la droite
                'Var7 = 1
                'Var5 = Length
                'Var8 = 0
                'Do While PosX > 0
                '    PosX = PosX - deltaX3
                '    If PosX < 0 Then
                '        deltaX3 = PosX
                '        PosX = 0
                '    End If
                '    Var1 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX) - Lambda(j, 3, Files + 1)) / Ksi(j, 3, Files + 1)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * Ksi(j, 3, Files + 1) * PosX)            'fs
                '    Var3 = System.Math.Exp(-0.5 * ((PosX - (Length - CDbl(TextBox1.Text))) / CDbl(TextBox2.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(TextBox2.Text))
                '    Var6 = System.Math.Exp(-0.5 * (((PosX - deltaX3 / 2) - (Length - CDbl(TextBox1.Text))) / CDbl(TextBox2.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(TextBox2.Text))
                '    Var2 = System.Math.Exp(-0.5 * ((Var5 - (Length - CDbl(TextBox1.Text))) / CDbl(TextBox2.Text)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * CDbl(TextBox2.Text))
                '    Var2 = (Var2 + Var6) * deltaX3 / 4         'méthode des trapèzes
                '    Var3 = (Var3 + Var6) * deltaX3 / 4 + Var2
                '    Var9 = Var1 * (Var8 + Var3) 'partie de Pf sans l'intégration de xi+1

                '    Var1 = System.Math.Exp(-0.5 * ((System.Math.Log(PosX - deltaX3 / 2) - Lambda(j, 3, Files + 1)) / Ksi(j, 3, Files + 1)) ^ 2) / ((2 * System.Math.PI) ^ 0.5 * Ksi(j, 3, Files + 1) * (PosX - deltaCL2 / 2))    'fs
                '    Var6 = Var1 * (Var8 + Var2) 'partie de Pf sans l'intégration de (xi+1 + xi)/2
                '    If PosX <> deltaX3 Then        'formule de simpson
                '        Var7 = (PosX - Var5) * (Var9 + 4 * Var6 + Var4) / 6
                '        LAen(j, 1) = LAen(j, 1) + Var7
                '    End If
                '    Var8 = Var8 + Var3
                '    Var4 = Var9
                '    Var5 = PosX
                'Loop
                'LAen(j, 1) = LAen(j, 1) / Var8  'correction due au troncage de l'épaisseur du mur
            Next j

            nfic = CShort(FreeFile())
            FileOpen(CInt(nfic), OT, OpenMode.Output)
            PrintLine(CInt(nfic), "Probabilité d'initiation de la corrosion due à la présence de ions chlorures")
            PrintLine(nfic, "temps, temps, Pf bord gauche")
            PrintLine(nfic, "années, jours,")
            For j = 1 To Nline + 2
                PrintLine(nfic, Lambda(j, 0, 1), ",", Lambda(j, 1, 1), ",", LAen(j, 0), ",")
            Next j
            FileClose(nfic)

            OT = "BCa_" & OutFile        'enregistrement de l'indice de fiabilité pour la corrosion par les chlorures
            ''''''''''''''''''''''''''''''''''
            Filtre = "txt files (BCa_*.txt)|BCa_*.txt|All files (*.*)|*.*"
            Index = CShort(1)
            Directoire = True
            Titre = "Indice de fiabilité (carb)"
            SaveDialog(OT, Canc, Filtre, Index, Directoire, Titre)
            If Canc = True Then GoTo f
            ''''''''''''''''''''''''''''''''''
            For j = 1 To Nline + 2  'calcul de l'indice de fiabilité bêta
                Var1 = -4           'beta
                Var2 = 4            'delta beta
                Var4 = 10
                Spos = False
                Spos_old = False
                Do While System.Math.Abs(Var2) > 0.001
                    If Var4 <> 10 Then Var1 = Var1 + Var2
                    If Var1 >= 0 Then
                        Var3 = 1 / (1 + 0.3275911 * Var1)
                        Var3 = 1 - (0.254829592 * Var3 - 0.284496736 * Var3 ^ 2 + 1.421413741 * Var3 ^ 3 - 1.453152027 * Var3 ^ 4 + 1.061405429 * Var3 ^ 5) * System.Math.Exp(-Var1 ^ 2)
                        Var3 = 0.5 * (1 - Var3)         'pf approximé
                    Else
                        Var3 = 1 / (1 + 0.3275911 * (-Var1))
                        Var3 = 1 - (0.254829592 * Var3 - 0.284496736 * Var3 ^ 2 + 1.421413741 * Var3 ^ 3 - 1.453152027 * Var3 ^ 4 + 1.061405429 * Var3 ^ 5) * System.Math.Exp(-Var1 ^ 2)
                        Var3 = 0.5 * (1 + Var3)         'pf approximé
                    End If
                    Var4 = LAen(j, 0) - Var3
                    If Var4 < 0 Then
                        Spos = False
                    Else
                        Spos = True
                    End If
                    If Spos <> Spos_old Then
                        Var2 = -Var2 / 10
                    End If
                    Spos_old = Spos
                    If LAen(j, 0) = 0 Then
                        Var1 = 20
                        Exit Do
                    End If
                    If Var1 < -20 Then
                        Var1 = -20
                        Exit Do
                    End If
                Loop
                KSen(j, 0) = Var1           'indice de fiabilité bêta
            Next j

            nfic = CShort(FreeFile())
            FileOpen(CInt(nfic), OT, OpenMode.Output)
            PrintLine(CInt(nfic), "Indice de fiabilité pour la corrosion due à la carbonatation")
            PrintLine(nfic, "temps, temps, Bêta,")
            PrintLine(nfic, "années, jours,")
            For j = 1 To Nline + 2
                PrintLine(nfic, Lambda(j, 0, 1), ",", Lambda(j, 1, 1), ",", KSen(j, 0), ",")
                LAen(j, 0) = 0
                LAen(j, 1) = 0
            Next j
            FileClose(nfic)
            Files = Files + 1
f:
        End If

            Beep()

    End Sub

    'activation de l'approche probabiliste pour l'enrobage
    Private Sub CheckBox8_CheckedChanged(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles CheckBox8.CheckedChanged
        If CheckBox8.Checked = True Then
            TextBox1.Text = 35
            TextBox1.Enabled = True
            TextBox2.Text = 23
            TextBox2.Enabled = True
        Else
            TextBox1.Text = ""
            TextBox1.Enabled = False
            TextBox2.Text = ""
            TextBox2.Enabled = False
        End If
    End Sub

End Class
