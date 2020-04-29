Module declarations
    Dim Picture_Name As String
    Dim Picture_Number As Short

    Dim frm As frmSplash

    'Demarrage du programme
    Public Sub main()
        Initialize_Data()
        frm = New frmSplash()
        frm.ShowDialog()
    End Sub

    'Auteurs du programme
    Public Sub Initialize_Data()
        'Store
        Picture_Name = "Guido & David"
        Picture_Number = CShort(0)
    End Sub

End Module
